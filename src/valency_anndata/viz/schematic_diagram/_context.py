from typing import Optional
from anndata import AnnData
from ._schematic import adata_structure_svg
from ._utils import _show_svg

class _SchematicDiagramContext:
    """
    Context manager that captures a snapshot of an AnnData object on entry,
    and automatically renders a diff schematic on exit.
    The snapshot only exists for the duration of the block.

    If a filename is provided, it is forwarded to the final rendered diagram
    when the context exits successfully.
    """

    def __init__(self, adata: AnnData, *, filename: Optional[str] = None):
        self.adata = adata
        self.filename = filename
        self.snapshot = None

    def __enter__(self):
        # Take a snapshot of the current state
        self.snapshot = self.adata.copy()
        return self.adata

    def __exit__(self, exc_type, exc, tb):
        # Render diff if no exception
        if exc_type is None:
            dwg = adata_structure_svg(self.adata, diff_from=self.snapshot)
            _show_svg(dwg, filename=self.filename)

        # Forget the snapshot: no global registration
        self.snapshot = None

        # Do not suppress exceptions
        return False
