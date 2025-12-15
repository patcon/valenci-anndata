import numpy as np
from anndata import AnnData
from typing import Literal, Optional

def impute(
    adata: AnnData,
    *,
    strategy: Literal["zero", "mean", "median"] = "mean",
    source_layer: Optional[str] = None,
    target_layer: Optional[str] = None,
    overwrite: bool = False,
) -> None:
    """
    Impute NaN values in an AnnData matrix and store the result in a layer.

    Parameters
    ----------
    adata
        AnnData object.
    strategy
        Imputation strategy. Currently supports:
        - "zero": replace NaNs with 0
        - "mean": column-wise mean
        - "median": column-wise median
    source_layer
        Layer to read from. If None, uses adata.X.
    target_layer
        Layer to write to. Defaults to "X_imputed_<strategy>".
    overwrite
        Whether to overwrite an existing target layer.
    """
    if target_layer is None:
        target_layer = f"X_imputed_{strategy}"

    if not overwrite and target_layer in adata.layers:
        return

    # Select source matrix
    if source_layer is None:
        X = adata.X
    else:
        X = adata.layers[source_layer]

    if X is None:
        raise ValueError("No source matrix available for imputation.")

    # Work on a copy
    X = np.asarray(X, dtype=float).copy()

    nan_mask = np.isnan(X)
    if not nan_mask.any():
        adata.layers[target_layer] = X
        return

    if strategy == "zero":
        X[nan_mask] = 0.0

    elif strategy in {"mean", "median"}:
        reducer = np.nanmean if strategy == "mean" else np.nanmedian
        col_stats = reducer(X, axis=0)

        # Replace NaNs column-wise
        rows, cols = np.where(nan_mask)
        X[rows, cols] = col_stats[cols]

    else:
        raise ValueError(f"Unknown imputation strategy: {strategy!r}")

    adata.layers[target_layer] = X
