from typing import Optional
from huggingface_hub import snapshot_download
import valency_anndata as val

def aufstehen(
    translate_to: Optional[str] = None,
):
    """
    Polis conversation of 33k+ Germans, run by political party Aufstehen.

    This is largest Polis conversation run as of now, in fall 2018.

    See: <https://compdemocracy.org/Case-studies/2018-germany-aufstehen/>

    The data is pulled from an archive at:
    <https://huggingface.co/datasets/patcon/polis-aufstehen-2018>

    Note
    ----

    This dataset has been augmented by merging `is-meta` and `is-seed` statement
    data (missing from the official CSV export) that were retreived from the
    Polis API. Specifically, `is-meta` is required in order to reproduce outputs
    of the Polis data pipeline.

    Attribution
    -----------

    Data was gathered using the Polis software (see:
    <https://compdemocracy.org/polis> and
    <https://github.com/compdemocracy/polis>) and is sub-licensed under CC BY
    4.0 with Attribution to The Computational Democracy Project. The data and
    more information about how the data was collected can be found at the
    following link: <https://pol.is/report/r6xd526vyjyjrj9navxrj>
    """
    export_dir = snapshot_download(
        repo_id="patcon/polis-aufstehen-2018",
        repo_type="dataset",
        # Suppress HF_TOKEN warning.
        token=False,
    )
    adata = val.datasets.polis.load(source=export_dir, translate_to=translate_to)

    return adata

