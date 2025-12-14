import valency_anndata as val

adata = val.datasets.polis("https://pol.is/report/r2dfw8eambusb8buvecjt") # small
# adata = val.datasets.polis("https://pol.is/report/r29kkytnipymd3exbynkd") # Chile
# adata = val.datasets.polis("https://pol.is/3ntrtcehas")

print(adata)
if adata.X is not None:
    print(adata.X.shape)
print(adata.X)
print(adata.uns["schema"])
print(adata.var)