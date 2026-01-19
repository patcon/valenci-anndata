# -*- coding: utf-8 -*-
"""
Polis 2.0 Statement Analysis Recipe

This module provides a recipe function for analyzing Polis statement data
using text embeddings, UMAP projection, and hierarchical clustering (EVOC).
"""

import numpy as np
import pandas as pd
from datetime import datetime
from typing import Optional, Dict, Any, TYPE_CHECKING
import logging

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from polismath_commentgraph.core import EmbeddingEngine, ClusteringEngine # type: ignore[import-not-found]


def recipe_polis2_statements(
    adata,
    # Text embedding parameters
    batch_size: int = 32,
    show_progress: bool = False,
    key_added_text_embed: str = "X_embed",
    
    # UMAP projection parameters
    umap_params: Optional[Dict[str, Any]] = None,
    key_added_umap: str = "X_umap_statements",
    
    # Clustering parameters
    num_layers: int = 4,
    evoc_params: Optional[Dict[str, Any]] = None,
    key_added_labels: str = "evoc_polis2",
    
    # Engine parameters
    embedding_engine: Optional["EmbeddingEngine"] = None,
    clustering_engine: Optional["ClusteringEngine"] = None,
    
    # Metadata storage
    store_metadata: bool = True,
    metadata_key: str = "polis2_metadata",
    
    inplace: bool = True,
) -> Optional[Any]:
    """
    Run the Polis 2.0 statement analysis pipeline.
    
    This recipe performs text embedding, dimensionality reduction via UMAP,
    and hierarchical clustering using EVOC on Polis statement data.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with statements in .var and voting matrix in .X.
        Must contain 'content' column in adata.var with statement text.
    batch_size : int, optional (default: 32)
        Batch size for text embedding processing.
    show_progress : bool, optional (default: False)
        Whether to show progress bar during embedding.
    key_added_text_embed : str, optional (default: "X_embed")
        Key in adata.varm where text embeddings will be stored.
    umap_params : dict, optional
        UMAP parameters. If None, uses default UMAPParameters.
        Can be generated from UMAPParameters.to_dict().
        Example: {"n_components": 2, "metric": "cosine", "n_neighbors": 15, "min_dist": 0.1}
    key_added_umap : str, optional (default: "X_umap_statements")
        Key in adata.varm where UMAP projection will be stored.
    num_layers : int, optional (default: 4)
        Number of hierarchical clustering layers to create.
    evoc_params : dict, optional
        EVOC clustering parameters. If None, uses default EVOCParameters.
        Can be generated from EVOCParameters.to_dict().
        Example: {"min_samples": 5, "min_cluster_size": 5}
    key_added_labels : str, optional (default: "evoc_polis2")
        Base key for cluster labels. Will store all layers in adata.varm[key_added_labels]
        as a 2D array (n_statements × num_layers), and the top-level (finest) clustering
        in adata.var[f"{key_added_labels}_top"] for convenience.
    embedding_engine : EmbeddingEngine, optional
        Custom embedding engine. If None, creates a new instance.
    clustering_engine : ClusteringEngine, optional
        Custom clustering engine. If None, creates a new instance.
    store_metadata : bool, optional (default: True)
        Whether to store processing metadata in adata.uns.
    metadata_key : str, optional (default: "polis2_metadata")
        Key in adata.uns where metadata will be stored.
    inplace : bool, optional (default: True)
        Whether to modify adata in place or return a copy.
    
    Returns
    -------
    AnnData or None
        Returns AnnData if inplace=False, otherwise returns None and modifies in place.
    
    Updates
    -------
    adata.varm[key_added_text_embed] : np.ndarray
        Text embeddings matrix (n_statements × embedding_dim)
    adata.varm[key_added_umap] : np.ndarray
        UMAP 2D projection (n_statements × 2)
    adata.varm[key_added_labels] : np.ndarray
        All cluster labels as 2D array (n_statements × num_layers)
    adata.var[f"{key_added_labels}_top"] : pd.Categorical
        Top-level (finest) cluster labels for convenience
    adata.uns[metadata_key] : dict
        Metadata about the processing (if store_metadata=True)
    
    Examples
    --------
    Basic usage:
    
    >>> import valency_anndata as val
    >>> adata = val.datasets.polis.load("https://pol.is/report/r3epuappndxdy7dwtvwpb")
    >>> val.tools.recipe_polis2_statements(adata)
    
    Advanced configuration:
    
    >>> from polismath_commentgraph.schemas.dynamo_models import UMAPParameters, EVOCParameters
    >>> umap_params = UMAPParameters(
    ...     n_components=2,
    ...     metric="cosine",
    ...     n_neighbors=15,
    ...     min_dist=0.1,
    ... )
    >>> evoc_params = EVOCParameters(
    ...     min_samples=5,
    ...     min_cluster_size=5,
    ... )
    >>> val.tools.recipe_polis2_statements(
    ...     adata,
    ...     umap_params=umap_params.to_dict(),
    ...     evoc_params=evoc_params.to_dict(),
    ...     num_layers=6,
    ... )
    """
    # Lazy import of heavy modules - only load when function is called
    try:
        from polismath_commentgraph.core import EmbeddingEngine, ClusteringEngine # type: ignore[import-not-found]
    except ImportError as e:
        raise ImportError(
            "The 'polismath_commentgraph' package is required for this function. "
            "Install it with: pip install valency-anndata[polis2]"
        ) from e
    
    # Handle inplace
    adata = adata if inplace else adata.copy()
    
    # Validate input
    if "content" not in adata.var.columns:
        raise ValueError(
            "adata.var must contain 'content' column with statement text. "
            f"Available columns: {adata.var.columns.tolist()}"
        )
    
    # Initialize engines
    if embedding_engine is None:
        logger.info("Initializing EmbeddingEngine...")
        embedding_engine = EmbeddingEngine()
    
    if clustering_engine is None:
        logger.info("Initializing ClusteringEngine...")
        # Apply EVOC parameters if provided
        if evoc_params is not None:
            clustering_engine = ClusteringEngine(**evoc_params)
        else:
            clustering_engine = ClusteringEngine()

    # For typing.
    assert embedding_engine is not None
    assert clustering_engine is not None
    
    # Extract texts
    texts = adata.var["content"].tolist()
    logger.info(f"Processing {len(texts)} statements")
    
    # Step 1: Generate text embeddings
    logger.info(f"Step 1/3: Generating text embeddings (batch_size={batch_size})...")
    embeddings = embedding_engine.embed_batch(
        texts=texts,
        batch_size=batch_size,
        show_progress=show_progress,
    )
    adata.varm[key_added_text_embed] = embeddings
    logger.info(f"  Embeddings shape: {embeddings.shape}")
    
    # Step 2: Project to 2D using UMAP
    logger.info("Step 2/3: Projecting to 2D with UMAP...")
    
    # Handle UMAP parameters
    if umap_params is not None:
        # Create a new clustering engine with custom UMAP parameters
        umap_clustering_engine = ClusteringEngine(
            umap_n_components=umap_params.get("n_components", 2),
            umap_n_neighbors=umap_params.get("n_neighbors", 15),
            umap_min_dist=umap_params.get("min_dist", 0.1),
            umap_metric=umap_params.get("metric", "cosine"),
        )
        umap_projection = umap_clustering_engine.project_to_2d(
            embeddings=embeddings,
        )
    else:
        umap_projection = clustering_engine.project_to_2d(
            embeddings=embeddings,
        )
    
    adata.varm[key_added_umap] = umap_projection
    logger.info(f"  UMAP projection shape: {umap_projection.shape}")
    
    # Step 3: Create hierarchical clustering layers
    logger.info(f"Step 3/3: Creating {num_layers} clustering layers...")
    
    cluster_layers = clustering_engine.create_clustering_layers(
        embeddings=embeddings,
        num_layers=num_layers,
    )
    
    # Store all cluster labels in varm as a 2D array (reversed order: finest to coarsest)
    cluster_array = np.column_stack([layer for layer in reversed(cluster_layers)])
    adata.varm[key_added_labels] = cluster_array
    
    # Also store top-level (finest granularity) in var for convenience
    adata.var[f"{key_added_labels}_top"] = pd.Categorical(cluster_array[:, 0])
    
    # Log cluster counts for each layer
    for i in range(cluster_array.shape[1]):
        n_clusters = len(np.unique(cluster_array[:, i][cluster_array[:, i] >= 0]))
        logger.info(f"  Layer {i}: {n_clusters} clusters")
    
    # Store metadata
    if store_metadata:
        metadata = {
            "processed_date": datetime.now().isoformat(),
            "num_statements": len(texts),
            "embedding_dim": embeddings.shape[1],
            "num_layers": num_layers,
            "cluster_counts": [
                int(len(np.unique(layer[layer >= 0])))
                for layer in cluster_layers
            ],
            "parameters": {
                "batch_size": batch_size,
                "key_text_embed": key_added_text_embed,
                "key_umap": key_added_umap,
                "key_labels": key_added_labels,
                "umap_params": umap_params,
                "evoc_params": evoc_params,
            }
        }
        
        # Add source information if available
        if "source" in adata.uns:
            if "base_url" in adata.uns["source"]:
                metadata["base_url"] = adata.uns["source"]["base_url"]
            if "conversation_id" in adata.uns["source"]:
                metadata["conversation_id"] = adata.uns["source"]["conversation_id"]
        
        adata.uns[metadata_key] = metadata
        logger.info(f"Metadata stored in adata.uns['{metadata_key}']")
    
    logger.info("✓ Polis 2.0 recipe completed successfully!")
    
    if not inplace:
        return adata