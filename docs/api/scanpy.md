Since we take inspiration from the `scanpy` tool in our data structure and conventions, many of its methods can be used on our vote matrix data.

For convenience, we make all `scanpy`'s methods available like so:

```py hl_lines="3"
import valency_anndata as val

val.scanpy.logging.print_header()
```

This is functionally equivalent to:

```py hl_lines="3"
import scanpy as sc

sc.logging.print_header()
```

See [scanpy](https://scanpy.readthedocs.io/en/stable/api/index.html) for more methods you can experiment with via the `val.scanpy` namespace.