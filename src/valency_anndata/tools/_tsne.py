import scanpy as sc
from anndata import AnnData
from typing import overload

from scanpy._utils.random import _LegacyRandom

@overload
def tsne( # pyright: ignore[reportInconsistentOverload]
    adata: AnnData,
    *,
    n_pcs: int | None = None,
    use_rep: str | None = None,
    perplexity: float = 30,
    metric: str = "euclidean",
    early_exaggeration: float = 12,
    learning_rate: float = 1000,
    random_state: _LegacyRandom = 0,
    use_fast_tsne: bool = False,
    n_jobs: int | None = None,
    key_added: str | None = None,
    copy: bool = False,
) -> AnnData | None: ...

def tsne(adata, **kwargs) -> AnnData | None:
    sc.tl.tsne(adata, **kwargs)