import scanpy as sc
from anndata import AnnData
from typing import overload

import numpy as np
from typing import Literal
from scanpy._utils.random import _LegacyRandom
_InitPos = Literal["paga", "spectral", "random"]

@overload
def umap( # pyright: ignore[reportInconsistentOverload]
    adata: AnnData,
    *,
    min_dist: float = 0.5,
    spread: float = 1.0,
    n_components: int = 2,
    maxiter: int | None = None,
    alpha: float = 1.0,
    gamma: float = 1.0,
    negative_sample_rate: int = 5,
    init_pos: _InitPos | np.ndarray | None = "spectral",
    random_state: _LegacyRandom = 0,
    a: float | None = None,
    b: float | None = None,
    method: Literal["umap", "rapids"] = "umap",
    key_added: str | None = None,
    neighbors_key: str = "neighbors",
    copy: bool = False,
) -> AnnData | None: ...

def umap(adata, **kwargs) -> AnnData | None:
    sc.tl.umap(adata, **kwargs)