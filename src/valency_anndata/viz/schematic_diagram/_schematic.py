import math
import svgwrite
from anndata import AnnData
from ._svg_blocks import draw_grid_block
from ._utils import diff_text_style


# -------------------
# Layer stacking configuration
# -------------------
LAYER_X_OFFSET = -5  # horizontal shift per layer (depth effect)
LAYER_Y_OFFSET = 16   # vertical shift per layer
LAYER_LABEL_Y_SPACING = 4  # additional vertical space for layer labels
FONT_SIZE = 12

COLORS = {
    "X": "#4fba6f",
    "layers": "#4fba6f",
    "obs": "#efc41c",
    "var": "#2c96c0",
    "obsm": "#ef9021",
    "varm": "#194c61",
    "varp": "#965ba5",
    "obsp": "#f15c5a",
}

def draw_layer_rect(
    dwg, x0, y0, cell, layer_index, layer_name, n_rows, n_cols, color="#e67e22", status: str | None = None
):
    """
    Draw a single layer as a flat rectangle behind X, optionally coloring
    the label based on diff status.

    Parameters
    ----------
    dwg : svgwrite.Drawing
    x0, y0 : float
        Base X/Y of the main X block
    cell : float
        Size of a single cell
    layer_index : int
        Which layer (0 = closest to X, 1 = next below)
    layer_name : str
        Name of the layer
    n_rows, n_cols : int
        Rows and columns of the layer
    color : str
        Outline color
    status : str | None
        Optional diff status for the layer label ('added', 'removed', or None)
    """
    width = n_cols * cell
    height = n_rows * cell

    # Apply configurable shifts for depth
    x_shift = LAYER_X_OFFSET * (layer_index + 1)
    y_shift = LAYER_Y_OFFSET * (layer_index + 1)

    rect_x = x0 + x_shift
    rect_y = y0 + y_shift

    # Draw rectangle behind X
    dwg.add(
        dwg.rect(
            insert=(rect_x, rect_y),
            size=(width, height),
            fill="white",
            stroke=color,
            stroke_width=2,
        )
    )

    # Label at inside bottom edge
    label_y = rect_y + height - LAYER_LABEL_Y_SPACING
    label_style = diff_text_style(status)
    dwg.add(
        dwg.text(
            layer_name,
            insert=(rect_x + 5, label_y),
            font_size=FONT_SIZE,
            font_family="sans-serif",
            **label_style,
        )
    )

# ------------------------------------------------------------
# AnnData → SVG with diff
# ------------------------------------------------------------
def adata_structure_svg(adata: AnnData, diff_from: AnnData | None = None):
    """
    Render a schematic SVG of an AnnData object, optionally highlighting
    differences in .obs, .var, .layers, and .obsm relative to a snapshot.

    Parameters
    ----------
    adata : AnnData
        The AnnData object to visualize.
    diff_from : AnnData | None
        Optional snapshot to compare against for highlighting differences.

    Returns
    -------
    svgwrite.Drawing
        The generated SVG drawing.
    """
    cell = 18
    max_cells = 10
    pad = 40
    line_height = 14
    obs_key_spacing = 15  # horizontal spacing between rotated keys

    # -------------------
    # Determine obs/var keys and order
    # -------------------
    if diff_from is not None:
        obs_keys = list(diff_from.obs.keys())
        obs_keys += [k for k in adata.obs.keys() if k not in obs_keys]

        var_keys = list(diff_from.var.keys())
        var_keys += [k for k in adata.var.keys() if k not in var_keys]

        obsm_keys = list(diff_from.obsm.keys())
        obsm_keys += [k for k in adata.obsm.keys() if k not in obsm_keys]
    else:
        obs_keys = list(adata.obs.keys())
        var_keys = list(adata.var.keys())
        obsm_keys = list(adata.obsm.keys())

    # -------------------
    # Compute diff status for obs, var, layers
    # -------------------
    obs_status: dict[str, str] = {}
    var_status: dict[str, str] = {}
    layer_status: dict[str, str] = {}
    obsm_status: dict[str, str] = {}

    if diff_from is not None:
        obs_prev, obs_now = set(diff_from.obs.keys()), set(adata.obs.keys())
        var_prev, var_now = set(diff_from.var.keys()), set(adata.var.keys())
        layer_prev, layer_now = set(diff_from.layers.keys()), set(adata.layers.keys())
        obsm_prev, obsm_now = set(diff_from.obsm.keys()), set(adata.obsm.keys())

        for key in obs_now - obs_prev:
            obs_status[key] = "added"
        for key in obs_prev - obs_now:
            obs_status[key] = "removed"

        for key in var_now - var_prev:
            var_status[key] = "added"
        for key in var_prev - var_now:
            var_status[key] = "removed"

        for key in layer_now - layer_prev:
            layer_status[key] = "added"
        for key in layer_prev - layer_now:
            layer_status[key] = "removed"

        for key in obsm_now - obsm_prev:
            obsm_status[key] = "added"
        for key in obsm_prev - obsm_now:
            obsm_status[key] = "removed"

    # -------------------
    # Determine matrix sizes
    # -------------------
    obs_cells = min(max_cells, math.ceil(math.sqrt(adata.n_obs)))
    var_cells = min(max_cells, math.ceil(math.sqrt(adata.n_vars)))
    obsm_cells = min(max_cells, math.ceil(math.sqrt(next(iter(adata.obsm.values())).shape[0]) if adata.obsm else 0))

    X_width = var_cells * cell
    X_height = obs_cells * cell
    # BUG: Why do we need to subtract 30px here to get everything to line up?
    # TODO: Make the obsm width adaptive to label length (need not be related to cells)
    obsm_width = obsm_cells * cell # - 30 # this breaks wide .X data layouts
    obsm_height = X_height  # match X height

    # -------------------
    # Var block height
    # -------------------
    var_block_height = max(60, len(var_keys) * line_height)

    # -------------------
    # Obs block width
    # -------------------
    min_obs_width = 60
    needed_obs_width = len(obs_keys) * obs_key_spacing
    obs_width = max(min_obs_width, needed_obs_width)

    tilt_factor = 0.707  # sin/cos 45°
    last_key_extra = len(obs_keys[-1]) * (FONT_SIZE * 0.5) * tilt_factor if obs_keys else 0
    extra_canvas_padding = last_key_extra + 10

    x0 = pad + 120
    y0 = pad + var_block_height + 30

    canvas_width = x0 + obsm_width + 20 + X_width + 30 + obs_width + extra_canvas_padding
    num_layers = len(adata.layers)
    max_layer_shift = LAYER_Y_OFFSET * num_layers
    canvas_height = y0 + X_height + max_layer_shift + 50

    dwg = svgwrite.Drawing(size=(canvas_width, canvas_height), profile="full")

    # -------------------
    # Layers (stacked behind X, preserve order)
    # -------------------
    all_layer_names = []
    if diff_from is not None:
        for ln in diff_from.layers.keys():
            all_layer_names.append(ln)
        for ln in adata.layers.keys():
            if ln not in all_layer_names:
                all_layer_names.append(ln)
    else:
        all_layer_names = list(adata.layers.keys())

    for i, layer_name in reversed(list(enumerate(all_layer_names))):
        layer_data = adata.layers.get(layer_name, None)
        if layer_data is not None:
            layer_rows = min(max_cells, math.ceil(math.sqrt(layer_data.shape[0])))
            layer_cols = min(max_cells, math.ceil(math.sqrt(layer_data.shape[1])))
        else:
            removed_data = diff_from.layers.get(layer_name) if diff_from else None
            if removed_data is not None:
                layer_rows = min(max_cells, math.ceil(math.sqrt(removed_data.shape[0])))
                layer_cols = min(max_cells, math.ceil(math.sqrt(removed_data.shape[1])))
            else:
                layer_rows = layer_cols = max_cells

        draw_layer_rect(
            dwg,
            x0=x0 + obsm_width + 20,  # shift to X + layers area
            y0=y0,
            cell=cell,
            layer_index=i,
            layer_name=layer_name,
            n_rows=layer_rows,
            n_cols=layer_cols,
            color=COLORS["layers"],
            status=layer_status.get(layer_name),
        )

    # -------------------
    # obsm block (left of X, stacked like layers)
    # -------------------
    if obsm_keys:
        base_x = x0
        base_y = y0  # align first entry with X

        n_rows = obsm_cells
        n_cols = var_cells

        for i, key in reversed(list(enumerate(obsm_keys))):
            layer_index = i  # first entry = 0

            # Compute offsets: first entry aligned with base_y, subsequent entries shifted
            x_shift = LAYER_X_OFFSET * (layer_index - 1)
            y_shift = LAYER_Y_OFFSET * (layer_index - 1)
            entry_x = base_x + x_shift
            entry_y = base_y + y_shift

            draw_layer_rect(
                dwg,
                x0=entry_x,
                y0=entry_y,
                cell=cell,
                layer_index=0,  # already applied shift manually
                layer_name=key,
                n_rows=n_rows,
                n_cols=n_cols,
                color=COLORS["obsm"],
                status=obsm_status.get(key),
            )

        # Draw generic obsm label in the middle of the first entry
        first_entry_x = base_x
        first_entry_y = base_y
        label_width = n_cols * cell
        label_height = n_rows * cell

        # Only count current obsm entries for the label
        if diff_from is not None:
            current_obsm_keys = [k for k in obsm_keys if obsm_status.get(k) != "removed"]
        else:
            current_obsm_keys = obsm_keys

        lines = ["obsm", f"{adata.n_obs} @ {len(current_obsm_keys)}"]
        line_height = 14
        start_y = first_entry_y + label_height / 2 - (len(lines) - 1) * line_height / 2
        for i, line in enumerate(lines):
            dwg.add(
                dwg.text(
                    line,
                    insert=(first_entry_x + label_width / 2, start_y + i * line_height),
                    text_anchor="middle",
                    font_size=12,
                    font_family="sans-serif",
                )
            )

    # -------------------
    # X block
    # -------------------
    draw_grid_block(
        dwg,
        x=x0 + obsm_width + 20,
        y=y0,
        width=X_width,
        height=X_height,
        rows=obs_cells,
        cols=var_cells,
        label=f"X\n{adata.n_obs} x {adata.n_vars}",
        stroke=COLORS["X"],
    )

    # -------------------
    # Obs block
    # -------------------
    draw_grid_block(
        dwg,
        x=x0 + obsm_width + 20 + X_width + 30,
        y=y0,
        width=obs_width,
        height=X_height,
        rows=obs_cells,
        cols=1,
        label=f"obs\n{adata.n_obs} x {adata.obs.shape[1]}",
        stroke=COLORS["obs"],
    )

    # -------------------
    # Obs keys (rotated)
    # -------------------
    baseline_y = y0 - 7
    for i, key in enumerate(obs_keys):
        x = x0 + obsm_width + 20 + X_width + 30 + 10 + i * obs_key_spacing
        style = diff_text_style(obs_status.get(key))
        dwg.add(
            dwg.text(
                key,
                insert=(x, baseline_y),
                font_size=FONT_SIZE,
                font_family="sans-serif",
                text_anchor="start",
                transform=f"rotate(-45,{x},{baseline_y})",
                **style,
            )
        )

    # -------------------
    # Var block (top)
    # -------------------
    draw_grid_block(
        dwg,
        x=x0 + obsm_width + 20,
        y=y0 - var_block_height - 30,
        width=X_width,
        height=var_block_height,
        rows=1,
        cols=var_cells,
        label=f"var\n{adata.n_vars} x {adata.var.shape[1]}",
        stroke=COLORS["var"],
    )

    # -------------------
    # Var keys
    # -------------------
    for i, key in enumerate(var_keys):
        y = y0 - 27 - (len(var_keys) - i) * line_height + line_height / 2
        style = diff_text_style(var_status.get(key))
        dwg.add(
            dwg.text(
                key,
                insert=(x0 + obsm_width + 20 - 10, y),
                font_size=FONT_SIZE,
                font_family="sans-serif",
                text_anchor="end",
                **style,
            )
        )

    return dwg