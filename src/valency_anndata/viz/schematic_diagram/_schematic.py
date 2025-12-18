import math
import svgwrite
from anndata import AnnData
from ._svg_blocks import draw_grid_block
from ._utils import diff_text_style

# ------------------------------------------------------------
# AnnData → SVG with diff
# ------------------------------------------------------------
def adata_structure_svg(adata: AnnData, diff_from: AnnData | None = None):
    cell = 18
    max_cells = 10
    pad = 40
    line_height = 14
    obs_key_spacing = 15  # horizontal spacing between rotated keys

    if diff_from is not None:
        # Start with the original order from diff_from
        obs_keys = list(diff_from.obs.keys())
        # Add any new keys that appear in adata but not in diff_from
        obs_keys += [k for k in adata.obs.keys() if k not in obs_keys]

        # Similarly for var keys
        var_keys = list(diff_from.var.keys())
        var_keys += [k for k in adata.var.keys() if k not in var_keys]
    else:
        obs_keys = list(adata.obs.keys())
        var_keys = list(adata.var.keys())

    # -------------------
    # Determine diff sets
    # -------------------
    obs_status: dict[str, str] = {}
    var_status: dict[str, str] = {}

    if diff_from is not None:
        obs_prev = set(diff_from.obs.keys())
        obs_now = set(adata.obs.keys())
        for key in obs_now - obs_prev:
            obs_status[key] = "added"
        for key in obs_prev - obs_now:
            obs_status[key] = "removed"

        var_prev = set(diff_from.var.keys())
        var_now = set(adata.var.keys())
        for key in var_now - var_prev:
            var_status[key] = "added"
        for key in var_prev - var_now:
            var_status[key] = "removed"

    # -------------------
    # Determine matrix size
    # -------------------
    obs_cells = min(max_cells, math.ceil(math.sqrt(adata.n_obs)))
    var_cells = min(max_cells, math.ceil(math.sqrt(adata.n_vars)))

    X_width = var_cells * cell
    X_height = obs_cells * cell

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

    font_size = 12
    tilt_factor = 0.707  # sin/cos 45°
    last_key_extra = len(obs_keys[-1]) * (font_size * 0.5) * tilt_factor if obs_keys else 0
    extra_canvas_padding = last_key_extra + 10

    x0 = pad + 120
    y0 = pad + var_block_height + 30
    canvas_width = x0 + X_width + 30 + obs_width + extra_canvas_padding
    canvas_height = X_height + var_block_height + 150

    dwg = svgwrite.Drawing(
        size=(canvas_width, canvas_height),
        profile="full",
    )

    # -------------------
    # X block
    # -------------------
    draw_grid_block(
        dwg,
        x=x0,
        y=y0,
        width=X_width,
        height=X_height,
        rows=obs_cells,
        cols=var_cells,
        label=f"X\n{adata.n_obs} x {adata.n_vars}",
        stroke="#2ecc71",
    )

    # -------------------
    # Obs block (right)
    # -------------------
    draw_grid_block(
        dwg,
        x=x0 + X_width + 30,
        y=y0,
        width=obs_width,
        height=X_height,
        rows=obs_cells,
        cols=1,
        label=f"obs\n{adata.n_obs} x {adata.obs.shape[1]}",
        stroke="#3498db",
    )

    # -------------------
    # Obs keys (rotated 45° above obs block)
    # -------------------
    baseline_y = y0 - 7
    for i, key in enumerate(obs_keys):
        x = x0 + X_width + 30 + 10 + i * obs_key_spacing
        style = diff_text_style(obs_status.get(key))

        dwg.add(
            dwg.text(
                key,
                insert=(x, baseline_y),
                font_size=font_size,
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
        x=x0,
        y=y0 - var_block_height - 30,
        width=X_width,
        height=var_block_height,
        rows=1,
        cols=var_cells,
        label=f"var\n{adata.n_vars} x {adata.var.shape[1]}",
        stroke="#9b59b6",
    )

    # -------------------
    # Var keys (left of var block)
    # -------------------
    for i, key in enumerate(var_keys):
        y = y0 - 27 - (len(var_keys) - i) * line_height + line_height / 2
        style = diff_text_style(var_status.get(key))

        dwg.add(
            dwg.text(
                key,
                insert=(x0 - 10, y),
                font_size=12,
                font_family="sans-serif",
                text_anchor="end",
                **style,
            )
        )

    return dwg