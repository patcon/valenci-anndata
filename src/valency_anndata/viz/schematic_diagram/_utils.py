import tempfile
import webbrowser

from ._browser import get_default_browser_name

def _display_svg_in_notebook(svg_text: str) -> bool:
    """Return True if displayed successfully in a notebook frontend."""
    try:
        from IPython import get_ipython # type: ignore[reportMissingImports]
        from IPython.display import SVG, display # type: ignore[reportMissingImports]
    except ImportError:
        return False

    ip = get_ipython()
    if ip is None:
        return False

    display(SVG(svg_text))
    return True

def _show_svg(dwg):
    svg_text = dwg.tostring()
    if _display_svg_in_notebook(svg_text):
        return

    # Otherwise, assume script mode → open in default browser
    with tempfile.NamedTemporaryFile(suffix=".svg", delete=False, mode="w") as f:
        f.write(svg_text)
        temp_svg_path = f.name
    browser = get_default_browser_name()
    if browser:
        webbrowser.get(browser).open(f"file://{temp_svg_path}")
        print(f"Opened in your default browser ({browser})")
    else:
        webbrowser.open(f"file://{temp_svg_path}")
        print("Opened in your browser")

def diff_text_style(status: str | None) -> dict[str, str]:
    if status == "added":
        return {
            "fill": "green",
            "font_weight": "bold",
        }
    if status == "removed":
        return {
            "fill": "red",
            "font_weight": "bold",
        }
    return {
        "fill": "black",
        "font_weight": "normal",
    }
