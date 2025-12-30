# Installation

For now, installation is a bit inelegant:

```python
pip install \
    git+https://github.com/patcon/valency-anndata \
    git+https://github.com/patcon/polis-client \
    git+https://github.com/polis-community/red-dwarf@algo-registries
```

Once we make a full release and settle the versions of our dependant libraries, it will become a simple `pip install valency-anndata`!

At the moment, there are a lot of core dependencies (even without installing development packages). In the future, we will will have a more slim core if there is interest.