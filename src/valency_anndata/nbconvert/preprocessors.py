from nbconvert.preprocessors import Preprocessor
import nbformat
import re

CODE_ANNOTATION_RE = re.compile(r"\(\d+\)")
ORDERED_LIST_RE = re.compile(r"^\s*\d+\.\s+", re.MULTILINE)

class MkdocsAnnotationPreprocessor(Preprocessor):
    """
    For code cells with mkdocs-material annotation markers:
    - clone into a code-only cell (no outputs)
    - insert the markdown annotation
    - clone outputs into a second empty code cell
    """

    def preprocess(self, nb, resources):
        new_cells = []
        pending_code = None

        for cell in nb.cells:
            if (
                cell.cell_type == "code"
                and CODE_ANNOTATION_RE.search(cell.source)
                and cell.outputs
            ):
                pending_code = cell
                continue

            if (
                pending_code
                and cell.cell_type == "markdown"
                and ORDERED_LIST_RE.match(cell.source or "")
            ):
                # Clone 1: code only
                code_only = nbformat.v4.new_code_cell(
                    source=pending_code.source,
                    metadata=pending_code.metadata,
                )

                # Clone 2: outputs only
                output_only = nbformat.from_dict(pending_code)
                output_only.metadata.setdefault("transient", {})
                output_only.metadata["transient"]["remove_source"] = True

                # Prevents annotation from including the next cell
                # if it happens to be text/plain formatted with `indent` processor.
                # TODO: Open pull request for nbconvert to use raw fenced code block instead of indent.
                # See: https://github.com/jupyter/nbconvert/blob/216550b2aae4c329f4dab597a96ae7cac30de79a/share/templates/markdown/index.md.j2#L36-L38
                md_spacer_cell = nbformat.v4.new_markdown_cell(source="---")

                new_cells.extend([
                    code_only,
                    cell,
                    md_spacer_cell,
                    output_only,
                ])

                pending_code = None
                continue

            if pending_code:
                # No annotation found; emit original cell unchanged
                new_cells.append(pending_code)
                pending_code = None

            new_cells.append(cell)

        if pending_code:
            new_cells.append(pending_code)

        nb.cells = new_cells
        return nb, resources
