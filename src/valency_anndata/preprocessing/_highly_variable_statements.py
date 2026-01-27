from anndata import AnnData
import numpy as np
import pandas as pd

def highly_variable_statements(
    adata: AnnData,
    *,
    layer: str | None = None,
    n_bins: int | None = 1,
    min_disp: float | None = None,
    max_disp: float | None = None,
    min_cov: int | None = 2,
    max_cov: int | None = None,
    n_top_statements: int | None = None,
    subset: bool = False,
    inplace: bool = True,
    variance_mode: str = "overall",  # "overall", "valence", "engagement"
    bin_by: str = "coverage",        # "coverage", "p_engaged", "mean_valence", "mean_abs_valence"
):
    """
    Identify highly variable statements in a vote matrix (AnnData).

    Parameters
    ----------
    variance_mode : str
        Which variance to use:
        - "overall" : variance of raw votes
        - "valence" : variance of engaged votes only
        - "engagement" : variance of engagement (1 if ±1, 0 if pass)
    bin_by : str
        What to bin on for normalization: "coverage", "p_engaged", "mean_valence", "mean_abs_valence"
    n_bins : int | None
        Number of bins for normalization; <=1 or None disables binning
    """

    # ---- 0. select matrix ---------------------------------------------
    X = adata.layers[layer] if layer is not None else adata.X
    X = np.asarray(X)
    n_statements = X.shape[1]

    # ---- 1. coverage and engagement -----------------------------------
    coverage = np.sum(~np.isnan(X), axis=0)
    engaged = (~np.isnan(X)) & (X != 0)
    p_engaged = engaged.sum(axis=0) / np.maximum(coverage, 1)

    # average valence for engaged votes only
    mean_valence = np.full(X.shape[1], np.nan)
    for j in range(X.shape[1]):
        vals = X[engaged[:, j], j]
        if vals.size > 0:
            mean_valence[j] = np.mean(vals)

    # optional: absolute version
    mean_abs_valence = np.abs(mean_valence)

    # ---- 2. compute variances -----------------------------------------
    # overall variance
    var_overall = np.nanvar(X, axis=0, ddof=1)

    # engagement variance: 1 if engaged, 0 if pass
    X_eng = np.where(np.isnan(X), np.nan, np.where(X != 0, 1.0, 0.0))
    var_engagement = np.nanvar(X_eng, axis=0, ddof=1)

    # valence variance: only consider engaged votes
    X_val = np.where(X == 0, np.nan, X)
    var_valence = np.nanvar(X_val, axis=0, ddof=1)

    # choose variance based on mode
    if variance_mode == "overall":
        dispersions = var_overall
    elif variance_mode == "valence":
        dispersions = var_valence
    elif variance_mode == "engagement":
        dispersions = var_engagement
    else:
        raise ValueError(f"Unknown variance_mode: {variance_mode}")

    valid = coverage >= 2  # same as before

    # ---- 3. binning ---------------------------------------------------
    if n_bins is None or n_bins <= 1:
        bin_idx = np.zeros(n_statements, dtype=int)
    else:
        if bin_by == "coverage":
            bin_idx = pd.cut(coverage, bins=n_bins, labels=False)
        elif bin_by == "mean_valence":
            bin_idx = pd.cut(mean_valence, bins=n_bins, labels=False)
        elif bin_by == "mean_abs_valence":
            bin_idx = pd.cut(mean_abs_valence, bins=n_bins, labels=False)
        elif bin_by == "p_engaged":
            bin_idx = pd.cut(p_engaged, bins=n_bins, labels=False)
        else:
            raise ValueError(f"Unknown bin_by: {bin_by}")

    # ---- 4. normalize within bins ------------------------------------
    dispersions_norm = np.full(n_statements, np.nan)
    for b in np.unique(bin_idx[valid]):
        mask = (bin_idx == b) & valid
        if mask.sum() < 2:
            continue
        d = dispersions[mask]
        mu = d.mean()
        sd = d.std()
        if sd == 0 or not np.isfinite(sd):
            continue
        dispersions_norm[mask] = (d - mu) / sd

    # ---- 5. stats table ----------------------------------------------
    stats = pd.DataFrame(
        {
            "coverage": coverage,
            "mean_valence": mean_valence,
            "mean_abs_valence": mean_abs_valence,
            "p_engaged": p_engaged,
            "bin_idx": bin_idx,
            "var_overall": var_overall,
            "var_valence": var_valence,
            "var_engagement": var_engagement,
            "dispersions": dispersions,
            "dispersions_norm": dispersions_norm,
        },
        index=adata.var_names,
    )

    # ---- 6. selection -------------------------------------------------
    if n_top_statements is not None:
        # rank by normalized dispersion first, then raw
        order = np.lexsort(
            (-stats["dispersions"].values, -stats["dispersions_norm"].values)
        )
        hv = np.zeros(n_statements, dtype=bool)
        hv[order[:n_top_statements]] = True
    else:
        hv = valid.copy()
        if min_cov is not None:
            hv &= stats["coverage"].values >= min_cov
        if max_cov is not None:
            hv &= stats["coverage"].values <= max_cov
        if min_disp is not None:
            hv &= stats["dispersions_norm"].values >= min_disp
        if max_disp is not None:
            hv &= stats["dispersions_norm"].values <= max_disp

    stats["highly_variable"] = hv

    # ---- 7. output ----------------------------------------------------
    if not inplace:
        return stats

    for k in stats.columns:
        adata.var[k] = stats[k].values

    # store metadata in .uns
    adata.uns["highly_variable"] = {
        "variance_mode": variance_mode,
        "bin_by": bin_by,
        "n_bins": n_bins,
        "min_disp": min_disp,
        "max_disp": max_disp,
        "min_cov": min_cov,
        "max_cov": max_cov,
        "n_top_statements": n_top_statements,
        "subset": subset,
        "valid": valid,
        "statement_names": adata.var_names.tolist(),
    }

    if subset:
        adata._inplace_subset_var(hv)
