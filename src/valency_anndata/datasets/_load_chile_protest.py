from typing import Optional
import valency_anndata as val

def chile_protest(
    translate_to: Optional[str] = None,
):
    """
    Polis conversation of 2,700+ Chileans during the 2019 #ChileDesperto protests.

    It was run informally by a single citizen, with minimal support
    infrastructure, outreach strategy, or moderation process.

    See: <https://en.wikipedia.org/wiki/Social_Outburst_(Chile)>

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
    following link: <https://pol.is/report/r29kkytnipymd3exbynkd>
    """
    adata = val.datasets.polis.load("https://pol.is/report/r29kkytnipymd3exbynkd", translate_to=translate_to)

    return adata

