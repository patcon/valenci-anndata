import valency_anndata as val
from valency_anndata.debug_utils import make_fake_polis_adata, fake_recipe_polis


USE_REAL_DATA = True

if USE_REAL_DATA:
    # adata = val.datasets.polis.load("https://pol.is/report/r2dfw8eambusb8buvecjt") # small
    adata = val.datasets.polis.load("https://pol.is/report/r29kkytnipymd3exbynkd", translate_to=None) # Chile
    # adata = val.datasets.polis.load("https://pol.is/3ntrtcehas")
    # adata = val.datasets.polis.load("https://pol.is/89euwydkce") # australia solar
    # adata = val.datasets.polis.load("https://pol.is/2demo") # minimum wage
else:
    adata = make_fake_polis_adata()

adata_snap = adata.copy()

with val.viz.schematic_diagram(diff_from=adata):
    if USE_REAL_DATA:
        val.tools.recipe_polis(
            adata,
            key_added_pca="X_pca",
            key_added_kmeans="kmeans_polis",
        )
    else:
        fake_recipe_polis(
            adata,
            key_added_pca="X_pca",
            key_added_kmeans="kmeans_polis",
        )

val.viz.schematic_diagram(adata)

print(adata_snap)
print("=========")
print(adata)

adata_snap_layers = adata.copy()
if not USE_REAL_DATA:
    del adata.obs["n_votes"]
del adata.layers["X_masked"]
del adata.obsm["X_pca_masked_unscaled"]
val.viz.schematic_diagram(adata, diff_from=adata_snap_layers)

val.scanpy.pl.pca(adata, color="kmeans_polis")