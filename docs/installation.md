# Installation

For now, installation is a bit inelegant:

```python
pip install \
    git+https://github.com/patcon/valency-anndata \
    git+https://github.com/patcon/polis-client \ # (1)!
    git+https://github.com/polis-community/red-dwarf@algo-registries # (2)!
```

1. [`polis-client`](https://github.com/patcon/polis-client) is an HTTP client library used to load data from Polis instance APIs.
2. [`red-dwarf`](https://github.com/polis-community/red-dwarf) is a utility library for functions related to the Polis data processing pipeline. General code in `valency-anndata` (code not reliant on `anndata`) will likely be upstreamed into `red-dwarf`.

Once we make a full release and settle the versions of our dependant libraries, installation will be simply `pip install valency-anndata`!

At the moment, there are a lot of core dependencies (even without installing development packages). In the future, we will have a more slim core if there is interest.