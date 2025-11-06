# Migration Guide: Old Code → New crested_utils

This guide shows how to migrate from your original scattered functions to the new unified `crested_utils` module.

## Table of Contents
1. [Basic Prediction](#basic-prediction)
2. [Cell Type Alignment](#cell-type-alignment)
3. [Correlation Calculation](#correlation-calculation)
4. [Complete Cross-Species Pipeline](#complete-cross-species-pipeline)

---

## Basic Prediction

### OLD CODE (from species_pred.py)
```python
from predict_chunked import predict_chunked

# Manual import and setup
adata_with_pred_species = predict_chunked(
    adata_species, 
    model, 
    chunk_size=1000,
    disk_based_saving=False
)
```

### NEW CODE
```python
from crested_utils import predict_regions

# Cleaner interface with better defaults
adata_with_pred_species = predict_regions(
    model=model,
    regions=adata_species,
    chunk_size=1000,
    verbose=True
)
```

**Benefits:**
- Clearer parameter names (`regions` instead of positional `adata`)
- Better documentation
- Consistent return value

---

## Cell Type Alignment

### OLD CODE (from species_pred.py)
```python
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
```

### NEW CODE
```python
from crested_utils import align_adata_cell_types

# All the above in 4 lines!
adata_species_aligned = align_adata_cell_types(
    adata_species,
    target_cell_types=celltypes_target,
    fill_missing=True,
    verbose=True
)
```

**Benefits:**
- **30+ lines → 4 lines**
- Handles all obsm and layers automatically
- Validated and tested
- Clear error messages
- Type-safe

---

## Combined: Prediction + Alignment

### OLD CODE
```python
# Two separate steps with manual code
adata_with_pred = predict_chunked(adata_species, model, chunk_size=2000, disk_based_saving=False)

# Then 30+ lines of manual alignment code...
celltypes_target = list(adata.obs_names)
present = [ct for ct in celltypes_target if ct in adata_species.obs_names]
# ... (all the alignment code)
```

### NEW CODE
```python
from crested_utils import predict_regions

# Single function does both!
adata_with_pred = predict_regions(
    model=model,
    regions=adata_species,
    target_cell_types=celltypes_target,
    align_cell_types=True,
    fill_missing_with_zeros=True,
    chunk_size=2000,
    verbose=True
)
```

**Benefits:**
- Single function call
- Alignment happens automatically
- No intermediate manual steps
- Preserved layers and metadata

---

## Correlation Calculation

### OLD CODE (from species_pred.py + rowwise_corr.py)
```python
import pandas as pd
from rowwise_corr import rowwise_corr
import matplotlib.pyplot as plt
import seaborn as sns

# Manual DataFrame creation
pred_data = pd.DataFrame(
    adata_with_pred_species.layers['pred_model'],
    index=adata_with_pred_species.obs_names
)
true_data = pd.DataFrame(
    adata_with_pred_species.X,
    index=adata_with_pred_species.obs_names
)

# Calculate correlation
corr = rowwise_corr(true_data.loc[common_celltypes], pred_data.loc[common_celltypes])

# Manual plotting
plt.figure(figsize=(20,20))
sns.heatmap(corr, cmap='Reds', vmin=0, vmax=1, square=True)
```

### NEW CODE
```python
from crested_utils import compare_predictions
import matplotlib.pyplot as plt
import seaborn as sns

# Single function for comparison
corr_matrix, self_corr = compare_predictions(
    adata_with_pred_species,
    prediction_layer='predicted',
    cell_types=common_celltypes,
    method='pearson'
)

# Get statistics
print(f"Mean self-correlation: {self_corr['correlation'].mean():.3f}")
print("\nPer cell type:")
print(self_corr.sort_values('correlation', ascending=False))

# Plotting
plt.figure(figsize=(12, 10))
sns.heatmap(corr_matrix, cmap='RdBu_r', center=0, vmin=-1, vmax=1, square=True)
plt.title('Predicted vs True Correlation')
plt.savefig('correlation.pdf')
```

**Benefits:**
- Built-in diagonal extraction
- Statistics included
- Cleaner interface
- Better default visualization

---

## Complete Cross-Species Pipeline

### OLD CODE (from species_pred.py)
```python
import anndata as ad
import crested
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from predict_chunked import predict_chunked
from rowwise_corr import rowwise_corr

species_common = ['chimpanzee']
genome_dir_map = {
    "chimpanzee": "panTro3",
}
base_genome = "/cluster/home/jjanssens/jjans/data/intestine/nhp_atlas/genomes/reference_"

for sp in species_common:
    print(f"=== {sp} ===")
    adata_path = f"../{sp}_crested/adata/{sp.capitalize()}_celltypes_filtered.h5ad"
    adata_species = ad.read_h5ad(adata_path)
    
    genome_folder = genome_dir_map[sp]
    genome_path = f"{base_genome}/{genome_folder}/fasta/genome.fa"
    genome = crested.Genome(genome_path)
    crested.register_genome(genome)
    
    common_celltypes = [x for x in human_celltypes if x in adata_species.obs_names]
    adata_species_filter = adata_species[common_celltypes,:].copy()
    
    # Optional filtering
    top_k = 500
    crested.pp.sort_and_filter_regions_on_specificity(
        adata_species_filter, top_k=top_k, method="proportion"
    )
    
    species_regions = list(set(adata_species_filter.var_names))
    adata_species = adata_species[common_celltypes, species_regions].copy()
    
    datamodule_species = crested.tl.data.AnnDataModule(
        adata_species,
        batch_size=64,
        genome=genome
    )
    
    # Prediction
    adata_with_pred_species = predict_chunked(
        adata_species, model, chunk_size=2000, disk_based_saving=False
    )
    
    # Manual alignment (30+ lines omitted for brevity)
    # ... alignment code ...
    
    # Manual correlation
    pred_data = pd.DataFrame(
        adata_with_pred_species.layers['pred_model'],
        index=adata_with_pred_species.obs_names
    )
    true_data = pd.DataFrame(
        adata_with_pred_species.X,
        index=adata_with_pred_species.obs_names
    )
    
    plt.figure(figsize=(20,20))
    sns.heatmap(
        rowwise_corr(true_data.loc[common_celltypes], pred_data.loc[common_celltypes]),
        cmap='Reds', vmin=0, vmax=1, square=True
    )
```

### NEW CODE
```python
import anndata as ad
import crested
import matplotlib.pyplot as plt
import seaborn as sns
from crested_utils import predict_regions, compare_predictions

species_list = ['chimpanzee']
genome_dir_map = {"chimpanzee": "panTro3"}
base_genome = "/cluster/home/jjanssens/jjans/data/intestine/nhp_atlas/genomes/reference_"

for sp in species_list:
    print(f"=== {sp} ===")
    
    # Load data
    adata_path = f"../{sp}_crested/adata/{sp.capitalize()}_celltypes_filtered.h5ad"
    adata_species = ad.read_h5ad(adata_path)
    
    # Setup genome
    genome_folder = genome_dir_map[sp]
    genome_path = f"{base_genome}/{genome_folder}/fasta/genome.fa"
    genome = crested.Genome(genome_path)
    crested.register_genome(genome)
    
    # Optional filtering
    common_celltypes = [x for x in human_celltypes if x in adata_species.obs_names]
    adata_species_filter = adata_species[common_celltypes, :].copy()
    
    crested.pp.sort_and_filter_regions_on_specificity(
        adata_species_filter, top_k=500, method="proportion"
    )
    
    # Predict with automatic alignment
    adata_pred = predict_regions(
        model=model,
        regions=adata_species_filter,
        target_cell_types=human_celltypes,
        align_cell_types=True,
        fill_missing_with_zeros=True,
        chunk_size=2000,
        verbose=True
    )
    
    # Compare predictions
    corr_matrix, self_corr = compare_predictions(
        adata_pred,
        cell_types=common_celltypes,
        method='pearson'
    )
    
    print(f"\nMean correlation: {self_corr['correlation'].mean():.3f}")
    
    # Plot
    plt.figure(figsize=(12, 10))
    sns.heatmap(
        corr_matrix,
        cmap='RdBu_r',
        center=0,
        vmin=-1,
        vmax=1,
        square=True
    )
    plt.title(f'{sp.capitalize()} - Predicted vs True')
    plt.savefig(f'{sp}_correlation.pdf')
```

**Benefits:**
- **~100 lines → ~40 lines** (60% reduction)
- Clearer logic flow
- No manual alignment code
- Built-in statistics
- Better error handling
- More maintainable

---

## Quick Reference Table

| Old Function/Code | New Function | Location |
|------------------|--------------|----------|
| `predict_chunked()` | `predict_regions()` | `crested_utils` |
| Manual alignment (30+ lines) | `align_adata_cell_types()` | `crested_utils` |
| `rowwise_corr()` | `rowwise_correlation()` | `crested_utils` |
| Manual DataFrame + correlation | `compare_predictions()` | `crested_utils` |
| `resize_region()` | `resize_region()` | `crested_utils` (unchanged) |

---

## Common Patterns

### Pattern 1: Quick Prediction
```python
# OLD
from predict_chunked import predict_chunked
result = predict_chunked(adata, model, chunk_size=1000, disk_based_saving=False)

# NEW
from crested_utils import predict_regions
result = predict_regions(model, adata, chunk_size=1000)
```

### Pattern 2: Cross-Species
```python
# OLD: 2 steps with manual code
result = predict_chunked(adata, model, ...)
# + 30 lines of alignment

# NEW: 1 step
result = predict_regions(
    model, adata, 
    target_cell_types=reference_celltypes,
    align_cell_types=True
)
```

### Pattern 3: Evaluation
```python
# OLD: Manual DataFrame creation
pred_df = pd.DataFrame(adata.layers['pred_model'], index=adata.obs_names)
true_df = pd.DataFrame(adata.X, index=adata.obs_names)
corr = rowwise_corr(true_df, pred_df)

# NEW: Built-in comparison
corr_matrix, self_corr = compare_predictions(adata)
print(f"Mean: {self_corr['correlation'].mean():.3f}")
```

---

## Testing Your Migration

After migrating, verify your results:

```python
# Test 1: Shape consistency
assert adata_new.shape == adata_old.shape

# Test 2: Prediction values close
import numpy as np
np.testing.assert_allclose(
    adata_new.layers['predicted'],
    adata_old.layers['pred_model'],
    rtol=1e-5
)

# Test 3: Cell type order
assert list(adata_new.obs_names) == list(adata_old.obs_names)

# Test 4: Correlation values similar
old_corr = old_method(...)
new_corr, _ = compare_predictions(...)
np.testing.assert_allclose(old_corr, new_corr, rtol=1e-5)
```

---

## Need Help?

1. Check `README.md` for full documentation
2. See `example_usage.py` for complete examples
3. Old function names still work (aliased for backward compatibility)
4. Start with one function at a time, don't migrate everything at once

---

## Summary

**Key Improvements:**
- ✅ 60% less code
- ✅ Automatic cell type alignment
- ✅ Built-in validation
- ✅ Better error messages
- ✅ Type hints for IDE support
- ✅ Comprehensive documentation
- ✅ Progress bars and logging
- ✅ Memory optimization

**Migration Effort:**
- Easy: Just replace function calls
- Most code: < 30 minutes
- Complex pipelines: < 2 hours
- Backward compatible: Old names still work
