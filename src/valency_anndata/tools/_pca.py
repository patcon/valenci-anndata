import scanpy as sc
from anndata import AnnData
from typing import overload
import numpy as np

from scanpy.preprocessing._pca import (
    SvdSolver,
    CSBase,
    _empty,
)
from scanpy._utils.random import _LegacyRandom
from scanpy._utils import Empty
from numpy.typing import DTypeLike, NDArray

@overload
def pca( # pyright: ignore[reportInconsistentOverload]
    data: AnnData | np.ndarray | CSBase,
    *,
    n_comps: int | None = None,
    layer: str | None = None,
    zero_center: bool = True,
    svd_solver: SvdSolver | None = None,
    chunked: bool = False,
    chunk_size: int | None = None,
    random_state: _LegacyRandom = 0,
    return_info: bool = False,
    mask_var: NDArray[np.bool_] | str | None | Empty = _empty,
    use_highly_variable: bool | None = None,
    dtype: DTypeLike = "float32",
    key_added: str | None = None,
    copy: bool = False,
) -> AnnData | np.ndarray | CSBase | None: ...

def pca(data, **kwargs) -> AnnData | np.ndarray | CSBase | None:
    sc.pp.pca(data, **kwargs)