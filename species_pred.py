import anndata as ad
import crested

species_common = ["macaque", "chimpanzee", "marmoset", "gorilla", "bonobo"]
species_common = ['chimpanzee']
# Map common name -> genome folder name
genome_dir_map = {
    "macaque": "Mmul10",
    "chimpanzee": "panTro3",
    "marmoset": "calJac1_mito",
    "gorilla": "gorGor4",
    "bonobo": "panPan1",
}

base_genome = "/cluster/home/jjanssens/jjans/data/intestine/nhp_atlas/genomes/reference_"

for sp in species_common:
    print(f"=== {sp} ===")
    # adata path: use common name for folder; capitalized for filename’s species token
    adata_path = f"../{sp}_crested/adata/{sp.capitalize()}_celltypes_filtered.h5ad"
    adata_species = ad.read_h5ad(adata_path)

    genome_folder = genome_dir_map[sp]
    genome_path = f"{base_genome}/{genome_folder}/fasta/genome.fa"

    genome = crested.Genome(genome_path)



    common_celltypes = [x for x in human_celltypes if x in adata_species.obs_names]
    adata_species_filter = adata_species[common_celltypes,:].copy()
    
    top_k = 500  # Here we take the top 2k most specific regions, but doing this on top 500 or top 1k will give similar results and will be faster to calculate
    crested.pp.sort_and_filter_regions_on_specificity(
        adata_species_filter, top_k=top_k, method="proportion"
    )
    
    
    species_regions = list(set(adata_species_filter.var_names))
    adata_species = adata_species[common_celltypes,species_regions].copy()
    
    datamodule_species = crested.tl.data.AnnDataModule(
            adata_species,
            batch_size=64,  # lower this if you encounter OOM errors
            genome=genome
        )
    
    
    crested.register_genome(
        genome
    )  # Register the genome so that it can be used by the package
    

    #adata_with_pred_species = predict_chunked(adata_species, model, chunk_size=1000,disk_based_saving=False)
import numpy as np
import pandas as pd
from anndata import AnnData

celltypes_target = list(adata.obs_names)  # desired final order

# Present vs missing in target order
present = [ct for ct in celltypes_target if ct in adata_species.obs_names]
missing = [ct for ct in celltypes_target if ct not in adata_species.obs_names]

# Subset to present
adata_present = adata_species[present, :].copy()

# --- Keep X dense float32 ---
X_present = np.asarray(adata_present.X)
if X_present.dtype != np.float32:
    X_present = X_present.astype(np.float32, copy=False)

# Build dense zeros for missing
n_missing = len(missing)
if n_missing:
    X_missing = np.zeros((n_missing, adata_species.n_vars), dtype=X_present.dtype)
    X_new = np.vstack([X_present, X_missing])
else:
    X_new = X_present

# Build obs in present+missing order (same as X_new)
obs_temp = pd.concat(
    [
        adata_present.obs,
        pd.DataFrame(index=pd.Index(missing, name=adata_species.obs_names.name)),
    ],
    axis=0,
)

# Copy var as-is
var_new = adata_species.var.copy()

# Create temp AnnData (present+missing order)
adata_temp = AnnData(X=X_new, obs=obs_temp, var=var_new)

# Preserve obsm['weights'] if available on the present set
if "weights" in adata_present.obsm:
    W = np.asarray(adata_present.obsm["weights"])
    if W.ndim == 1:
        W = W.reshape(-1, 1)
    W_new = np.vstack([W, np.zeros((n_missing, W.shape[1]), dtype=W.dtype)]) if n_missing else W
    adata_temp.obsm["weights"] = W_new

# Finally, reorder to exact human cell type order (affects X, obs, obsm consistently)
adata_species = adata_temp[celltypes_target, :].copy()

# Sanity: confirm dtype/shape
assert adata_species.X.dtype == np.float32
print(adata_species.X.shape, adata_species.X.dtype)

adata_with_pred_species = predict_chunked(adata_species, model, chunk_size=2000,disk_based_saving=False)

pred_data = pd.DataFrame(adata_with_pred_species.layers['pred_model'],index=adata_with_pred_species.obs_names)
true_data = pd.DataFrame(adata_with_pred_species.X,index=adata_with_pred_species.obs_names)


plt.figure(figsize=(20,20))
sns.heatmap(rowwise_corr(true_data.loc[common_celltypes],pred_data.loc[common_celltypes]),cmap='Reds',vmin=0,vmax=1,square=True)