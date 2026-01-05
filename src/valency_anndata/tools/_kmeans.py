import numpy as np
import pandas as pd
from anndata import AnnData
from reddwarf.sklearn.cluster import BestPolisKMeans
from typing import Literal, Optional, Tuple
from scanpy.get import _check_mask
from scanpy.tools._utils import _choose_representation
from numpy.typing import NDArray


def kmeans(
    adata: AnnData,
    use_rep: Optional[str] = None,
    n_pcs: Optional[int] = None,
    k_bounds: Optional[Tuple[int, int]] = None,
    init: Literal["kmeans++", "random", "polis"] = "kmeans++",
    init_centers: Optional[np.ndarray] = None,
    random_state: Optional[int] = None,
    mask_obs: NDArray[np.bool_] | str | None = None,
    key_added: str = "kmeans",
    inplace: bool = True,
) -> AnnData | None:
    """
    Apply BestPolisKMeans clustering to an AnnData object.

    Parameters
    ----------
    adata :
        Input data. Must have `.X` as a numpy array.
    use_rep
        Representation to use for clustering. If ``None``, use ``'X_pca'`` if
        present in ``adata.obsm``, otherwise fall back to ``adata.X``.
    n_pcs
        Number of dimensions to use from the selected representation. If given,
        only the first ``n_pcs`` columns are used.
    k_bounds :
        Minimum and maximum number of clusters to try. Defaults to [2, 5].
    init :
        Initialization method for KMeans. Defaults to 'polis'.
    init_centers :
        Initial cluster centers to use.
    random_state :
        Random seed for reproducibility.
    mask_obs :
        Restrict clustering to a certain set of observations. The mask is
        specified as a boolean array or a string referring to an array in
        [anndata.AnnData.obs][].
    key_added :
        Name of the column to store cluster labels in `adata.obs`.
    inplace :
        If True, modify `adata` in place and return None.
        If False, return a copy with the clustering added.

    Returns
    -------
    AnnData or None
        Returns a copy if `inplace=False`, otherwise modifies in place.
    """
    X = _choose_representation(adata, use_rep=use_rep, n_pcs=n_pcs)

    if not isinstance(X, np.ndarray):
        raise ValueError("Selected representation must be a numpy array.")

    if k_bounds is None:
        k_bounds_list = [2, 5]
    else:
        k_bounds_list = list(k_bounds)

    mask = _check_mask(adata, mask_obs, "obs")
    if mask is None:
        X_cluster = X
    else:
        X_cluster = X[mask]
        if X_cluster.shape[0] == 0:
            raise ValueError("mask_obs excludes all observations.")

    best_kmeans = BestPolisKMeans(
        k_bounds=k_bounds_list,
        init=init,
        init_centers=init_centers,
        random_state=random_state,
    )
    best_kmeans.fit(X_cluster)

    if not best_kmeans.best_estimator_:
        raise RuntimeError("BestPolisKMeans did not find a valid estimator.")

    raw_labels = best_kmeans.best_estimator_.labels_

    if mask is None:
        full_labels = raw_labels
    else:
        # dtype=object keeps labels from casting to float.
        full_labels = np.full(adata.n_obs, np.nan, dtype=object)
        full_labels[mask] = raw_labels

    labels = pd.Categorical(full_labels)

    def _write_kmeans_result(adata_out: AnnData) -> None:
        adata_out.obs[key_added] = labels

        kmeans_params = dict(
            k_bounds=k_bounds_list,
            best_k=best_kmeans.best_k_,
            best_score=best_kmeans.best_score_,
            init=init,
            random_state=random_state,
            use_rep=use_rep,
            n_pcs=n_pcs,
        )

        adata_out.uns[key_added] = {}
        adata_out.uns[key_added]["params"] = kmeans_params

    if inplace:
        _write_kmeans_result(adata)
        return None
    else:
        adata_copy = adata.copy()
        _write_kmeans_result(adata_copy)
        return adata_copy