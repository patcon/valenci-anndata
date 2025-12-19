from anndata import AnnData
import numpy as np
import pandas as pd

def make_fake_polis_adata(
    n_participants: int = 400,
    n_statements: int = 60,
) -> AnnData:
    """
    Create a fake Polis-like AnnData object with:
      - X: vote matrix (participants × statements)
      - obs: participant metadata
      - var: statement metadata
    """

    rng = np.random.default_rng(42)

    # -------------------
    # Fake vote matrix
    # Values roughly match Polis: -1 disagree, 0 pass, 1 agree
    # -------------------
    X = rng.choice(
        [-1.0, 0.0, 1.0],
        size=(n_participants, n_statements),
        p=[0.35, 0.30, 0.35],
    )

    # -------------------
    # obs: participants
    # -------------------
    obs = pd.DataFrame(
        {
            "participant_id": [f"p{i:04d}" for i in range(n_participants)],
            "n_votes": (X != 0).sum(axis=1),
        }
    ).set_index("participant_id")

    # -------------------
    # var: statements
    # -------------------
    var = pd.DataFrame(
        {
            "statement_id": [f"s{i:03d}" for i in range(n_statements)],
            "text_length": rng.integers(20, 200, size=n_statements),
        }
    ).set_index("statement_id")

    return AnnData(X=X, obs=obs, var=var)

def fake_recipe_polis(
    adata: AnnData,
    *,
    key_added_pca: str = "X_pca",
    key_added_kmeans: str = "kmeans_polis",
    n_pca: int = 2,
    n_clusters: int = 5,
):
    rng = np.random.default_rng(0)
    n_obs = adata.n_obs
    n_vars = adata.n_vars

    # -------------------
    # Fake PCA (obsm)
    # -------------------
    # Add the requested obsm entries
    adata.obsm["X_pca_masked_unscaled"] = np.zeros((n_obs, n_pca))
    adata.obsm["X_pca_masked_scaled"] = np.zeros((n_obs, n_pca))
    adata.obsm[key_added_pca] = rng.normal(size=(n_obs, n_pca))

    # -------------------
    # Fake varm entries
    # -------------------
    adata.varm["X_pca_masked_unscaled"] = np.zeros((n_vars, n_pca))

    # -------------------
    # Fake KMeans labels (obs)
    # -------------------
    labels = rng.integers(0, n_clusters, size=n_obs)
    adata.obs[key_added_kmeans] = pd.Categorical(labels)

    # -------------------
    # Fake layers
    # -------------------
    X_masked = adata.X.copy()
    mask = rng.random(size=X_masked.shape) < 0.2
    X_masked[mask] = 0
    adata.layers["X_masked"] = X_masked

    X_imputed = X_masked.copy()
    col_means = np.where(X_masked != 0, X_masked, np.nan).mean(axis=0)
    inds = np.where(X_imputed == 0)
    X_imputed[inds] = np.take(col_means, inds[1])
    adata.layers["X_imputed_mean"] = X_imputed

    # -------------------
    # Provenance marker (uns)
    # -------------------
    adata.uns.setdefault("polis", {})
    adata.uns["polis"]["fake"] = True