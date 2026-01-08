import ipywidgets as widgets
from IPython.display import display
from jscatter import Scatter as JScatter
from typing import Sequence, Tuple, Iterable
import pandas as pd
from anndata import AnnData


def obsm_to_df(
    adata,
    projections: Sequence[Tuple[str, str]],
    n_dims: int = 2,
    obs_cols: Iterable[str] | None = None,
) -> pd.DataFrame:
    dfs = []

    for key, prefix in projections:
        X = adata.obsm[key]
        cols = [f"{prefix}{i+1}" for i in range(n_dims)]
        dfs.append(
            pd.DataFrame(
                X[:, :n_dims],
                index=adata.obs_names,
                columns=pd.Index(cols),
            )
        )

    if obs_cols:
        dfs.append(adata.obs[list(obs_cols)].copy())

    return pd.concat(dfs, axis=1)

# TODO: Single representation case: no button
# TODO: kind="single" vs kind="multi"
def jscatter(
    adata: AnnData,
    use_reps: list[str] = [],
    color: str | None = None,
    height: int = 640,
    dark_mode: bool = True,
) -> JScatter:
    """
    Interactive Jupyter-Scatter view over one or more embeddings.

    Parameters
    ----------

    adata :
        An AnnData object with some projected representations stored in
        [`.obsm`][anndata.AnnData.obsm].
    use_reps :
        Projected representations of the data stored in [`.obsm`][anndata.AnnData.obsm].
    color :
        Key in [`.obs`][anndata.AnnData.obs] for coloring projections.
    height :
        Height of the scatter widget in cell output.
    dark_mode :
        Whether to set the plot background dark.

    Returns
    -------

    scatter : JScatter
        A populated [`Scatter`](https://jupyter-scatter.dev/api#scatter) instance.

    Examples
    --------

    ```py
    scatter = val.viz.jscatter(adata, use_reps=["X_pca_polis", "X_pacmap"], color="kmeans_polis")
    scatter.show()
    ```
    """
    background = "#1E1E20" if dark_mode else None

    # ---- prepare projections ----
    projections = [
        (key, key.removeprefix("X_").split("_")[0])
        for key in use_reps
    ]

    obs_cols = [color] if color is not None else None

    df = obsm_to_df(
        adata,
        projections=projections,
        obs_cols=obs_cols,
    )

    # ---- create scatter ----
    default_prefix = projections[0][1]

    scatter = JScatter(
        data=df,
        x=f"{default_prefix}1",
        y=f"{default_prefix}2",
        height=height,
    )

    if color is not None:
        scatter.color(by=color)

    scatter.background(background)

    # ---- projection toggle ----
    toggle = widgets.ToggleButtons(
        options=[
            (prefix.upper(), prefix)
            for _, prefix in projections
        ],
        value=default_prefix,
        description="Projection:",
    )

    def on_toggle(change):
        prefix = change["new"]
        scatter.xy(f"{prefix}1", f"{prefix}2")

    toggle.observe(on_toggle, names="value")

    display(toggle)

    return scatter
