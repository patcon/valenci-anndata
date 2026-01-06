import numpy as np
from reddwarf.sklearn.transformers import calculate_scaling_factors
from anndata import AnnData
import valency_anndata as val

def _zero_mask(
    adata: AnnData,
    *,
    key_added_var_mask: str = "zero_mask",
    key_added_layer: str = "X_masked",
    inplace: bool = True,
) -> AnnData | None:
    if not inplace:
        adata = adata.copy()

    if bool(adata.var["is_meta"].isna().any()):
        raise ValueError("Your statements are missing is_meta data. Cannot create required zero_mask. Force all values to be False to override.")

    adata.var[key_added_var_mask] = adata.var.eval("~is_meta and moderation_state > -1")

    mask = adata.var[key_added_var_mask].to_numpy()

    # Preconditions
    assert isinstance(adata.X, np.ndarray)

    X_masked = adata.X.copy()
    X_masked[:, ~mask] = 0
    adata.layers[key_added_layer] = X_masked

    if not inplace:
        return adata

def _sparsity_aware_scaling(
    adata: AnnData,
    use_rep: str = "X_pca_masked_unscaled",
    key_added: str = "X_pca_masked_scaled",
) -> AnnData | None:
    # NOTE: This calculates row sparsity based on zero'd out values
    # (meta & moderated out). This feels strange, but is consistent
    # with what Polis does.
    X_sparse = adata.layers["X_masked"]

    # Preconditions
    assert isinstance(X_sparse, np.ndarray)

    # 4. Scale PCA using sparsity data
    scaling_factors = calculate_scaling_factors(X_sparse)
    X_pca_unscaled = adata.obsm[use_rep]
    adata.obsm[key_added] = X_pca_unscaled / scaling_factors[:, None]

def _cluster_mask(
    adata: AnnData,
    participant_vote_threshold: int = 7,
    key_added_obs_mask: str = "cluster_mask",
) -> AnnData | None:
    # Preconditions
    assert isinstance(adata.X, np.ndarray)

    non_nan_counts = np.sum(~np.isnan(adata.X), axis=1)
    adata.obs[key_added_obs_mask] = non_nan_counts >= participant_vote_threshold

def recipe_polis(
    adata: AnnData,
    *,
    participant_vote_threshold: int = 7,
    key_added_pca: str = "X_pca_polis",
    key_added_kmeans: str = "kmeans_polis",
    inplace: bool = True,
):
    if not inplace:
        adata = adata.copy()

    # Preconditions
    assert isinstance(adata.X, np.ndarray)

    # 1. Mask statements with zeros
    _zero_mask(
        adata,
        key_added_var_mask="zero_mask",
        key_added_layer="X_masked",
    )

    # 2. Impute
    val.preprocessing.impute(
        adata,
        strategy="mean",
        source_layer="X_masked",
        target_layer="X_masked_imputed_mean",
    )

    # 3. PCA (unscaled)
    val.tools.pca(
        adata,
        layer="X_masked_imputed_mean",
        key_added="X_pca_masked_unscaled",
    )

    # 4. Scale PCA using sparsity data
    _sparsity_aware_scaling(
        adata,
        use_rep="X_pca_masked_unscaled",
        key_added=key_added_pca,
    )

    # Create cluster mask for threshold
    _cluster_mask(
        adata,
        participant_vote_threshold=participant_vote_threshold,
        key_added_obs_mask="cluster_mask",
    )

    # 5. KMeans clustering
    val.tools.kmeans(
        adata,
        use_rep=key_added_pca,
        # Force kmeans to only run on first two principle components.
        n_pcs=2,
        k_bounds=(2, 5),
        init="polis",
        mask_obs="cluster_mask",
        key_added=key_added_kmeans,
        inplace=inplace,
    )

    if not inplace:
        return adata