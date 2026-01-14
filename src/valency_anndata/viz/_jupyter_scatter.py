import ipywidgets as widgets
from IPython.display import display
from jscatter import Scatter as JScatter, compose
from typing import Optional, Sequence, Tuple, Iterable
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
    color: str | Iterable[str] | None = None,
    height: int = 640,
    dark_mode: bool = True,
    nrows: Optional[int] = None,
    ncols: Optional[int] = None,
    return_objs: bool = False,
) -> list[JScatter] | None:
    """
    Interactive Jupyter-Scatter view showing one or more embeddings. [[Lekschas _et al._, 2024](https://doi.org/10.21105/joss.07059)]

    A button is created for each projected representation, and clicking will
    animate points into that projection.

    Passing multiple color keys will display mulitple linked views.

    Parameters
    ----------

    adata :
        An AnnData object with some projected representations stored in
        [`.obsm`][anndata.AnnData.obsm].
    use_reps :
        One or more keys for projected representations of the data stored in
        [`.obsm`][anndata.AnnData.obsm].
    color :
        One or more keys in [`.obs`][anndata.AnnData.obs] for coloring each participant.
        Categorical values will use a discrete color map
        ([`okabeito`](https://cmap-docs.readthedocs.io/en/latest/catalog/qualitative/okabeito:okabeito/)),
        and anything else will use a continuous gradient
        ([`viridis`](https://cmap-docs.readthedocs.io/en/latest/catalog/sequential/bids:viridis/)).
    height :
        Pixel height of the scatter widget in output cell.
    dark_mode :
        Whether to set the plot background dark.
    nrows :
        Number of rows to display the scatter plots in.
    ncols :
        Number of columns to display the scatter plots in.
    return_objs :
        Whether to return the Scatter object(s).

    Returns
    -------

    scatters :
        A list of [`Scatter`](https://jupyter-scatter.dev/api#scatter) instances.

    Examples
    --------

    Plotting multiple representations in one view, colored with discrete categorical values.

    ```py
    val.viz.jscatter(
        adata,
        use_reps=["X_pca_polis", "X_localmap"],
        color="kmeans_polis",
    )
    ```

    <img src="../../assets/documentation-examples/viz--jscatter--single.png">

    Plotting mulitple `.obs` keys across multiple views, colored with continuous values.

    ```py
    val.viz.jscatter(
        adata,
        use_reps=["X_pca_polis", "X_pacmap"],
        color=["n_votes", "pct_agree", "pct_pass", "pct_disagree"],
        height=320,
    )
    ```

    <img src="../../assets/documentation-examples/viz--jscatter--multi.png">
    """
    background = "#1E1E20" if dark_mode else None

    # ---- prepare projections ----
    projections = [
        (key, key.removeprefix("X_").split("_")[0])
        for key in use_reps
    ]

    if color is None:
        colors = []
    elif isinstance(color, str):
        colors = [color]
    else:
        colors = list(color)

    obs_cols = colors if colors else None

    df = obsm_to_df(
        adata,
        projections=projections,
        obs_cols=obs_cols,
    )

    # ---- create scatter(s) ----
    default_prefix = projections[0][1]

    scatters = []

    for c in colors or [None]:
        scatter = JScatter(
            data=df,
            x=f"{default_prefix}1",
            y=f"{default_prefix}2",
            height=height,
            zoom_on_selection=True,
        )

        if c is not None:
            scatter.color(by=c)

        scatter.background(background)

        scatters.append(scatter)

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
        for s in scatters:
            s.xy(f"{prefix}1", f"{prefix}2")

    toggle.observe(on_toggle, names="value")

    grid = compose(
        list(zip(scatters, colors)),
        sync_view=True,
        sync_selection=True,
        sync_hover=True,
        cols=ncols,
        rows=nrows,
        row_height=height,
    )

    display(toggle)
    display(grid)

    return scatters if return_objs else None
