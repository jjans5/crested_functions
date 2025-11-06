# Quick Reference: crested_utils

## Installation

```bash
pip install crested anndata pandas numpy
pip install psutil  # Optional: for memory monitoring
```

## Basic Import

```python
from crested_utils import (
    predict_regions,           # Main prediction function
    align_adata_cell_types,    # Align cell types
    compare_predictions,       # Compare pred vs true
    rowwise_correlation,       # Row-wise correlation
    resize_region              # Resize regions
)
```

---

## Function Cheat Sheet

### 1. predict_regions() - Main Function

```python
# Basic prediction
adata_pred = predict_regions(model, regions=adata)

# With cell type alignment (for cross-species)
adata_pred = predict_regions(
    model=model,
    regions=adata_species,
    target_cell_types=['T_cell', 'B_cell', 'NK_cell', ...],
    align_cell_types=True,
    fill_missing_with_zeros=True,
    chunk_size=1000,
    batch_size=16,
    layer_name='predicted',
    verbose=True
)
```

**Parameters:**
- `model` - Your trained CREsted model
- `regions` - AnnData or list of region strings
- `target_cell_types` - Cell types for alignment (optional)
- `align_cell_types` - Align to target (default: True)
- `fill_missing_with_zeros` - Fill missing cell types (default: True)
- `chunk_size` - Regions per chunk (default: 1000)
- `layer_name` - Output layer name (default: 'predicted')

**Returns:** AnnData with predictions in `.layers[layer_name]`

---

### 2. align_adata_cell_types() - Cell Type Alignment

```python
# Add missing cell types as zeros
adata_aligned = align_adata_cell_types(
    adata,
    target_cell_types=['T_cell', 'B_cell', ...],
    fill_missing=True
)

# Or keep only common cell types
adata_aligned = align_adata_cell_types(
    adata,
    target_cell_types=['T_cell', 'B_cell', ...],
    fill_missing=False
)
```

**Parameters:**
- `adata` - Input AnnData
- `target_cell_types` - Desired cell types in order
- `fill_missing` - Fill (True) or drop (False) missing cell types

**Returns:** New AnnData with aligned cell types

---

### 3. compare_predictions() - Evaluate Predictions

```python
# Compare predictions vs true values
corr_matrix, self_corr = compare_predictions(
    adata,
    prediction_layer='predicted',
    true_layer=None,  # None = use .X
    cell_types=None,  # None = use all
    method='pearson'  # or 'spearman'
)

# View results
print(f"Mean correlation: {self_corr['correlation'].mean():.3f}")
print(self_corr.sort_values('correlation', ascending=False))
```

**Parameters:**
- `adata` - AnnData with predictions
- `prediction_layer` - Prediction layer name
- `true_layer` - True values layer (None = .X)
- `cell_types` - Cell types to compare (None = all)
- `method` - 'pearson' or 'spearman'

**Returns:**
- `corr_matrix` - Full correlation matrix (n x n)
- `self_corr` - DataFrame with diagonal correlations + cell types

---

### 4. rowwise_correlation() - Direct Correlation

```python
import pandas as pd

# Create DataFrames
pred_df = pd.DataFrame(adata.layers['predicted'], index=adata.obs_names)
true_df = pd.DataFrame(adata.X, index=adata.obs_names)

# Compute correlation
corr = rowwise_correlation(pred_df, true_df, method='pearson')
```

**Parameters:**
- `A` - First DataFrame (n_rows_A x n_columns)
- `B` - Second DataFrame (n_rows_B x n_columns)
- `method` - 'pearson' or 'spearman'

**Returns:** Correlation matrix (n_rows_A x n_rows_B)

---

### 5. resize_region() - Resize Genomic Regions

```python
# Resize to fixed length (centered on midpoint)
new_region = resize_region('chr1:1000-2000', new_length=500)

# Resize centered on specific position
new_region = resize_region('chr1:1000-2000', new_length=500, summit=1600)
```

**Parameters:**
- `region` - Region string ('chr:start-end')
- `new_length` - Desired length
- `summit` - Center position (None = midpoint)

**Returns:** Resized region string

---

## Common Workflows

### Workflow 1: Simple Prediction

```python
import crested
import anndata as ad
from crested_utils import predict_regions

# Load
model = crested.load_model('model.h5')
adata = ad.read_h5ad('data.h5ad')

# Predict
adata_pred = predict_regions(model, regions=adata)

# Access
predictions = adata_pred.layers['predicted']
```

---

### Workflow 2: Cross-Species Comparison (YOUR MAIN USE CASE)

```python
from crested_utils import predict_regions, compare_predictions
import matplotlib.pyplot as plt
import seaborn as sns

# Define reference cell types
human_celltypes = ['T_cell', 'B_cell', 'NK_cell', 'Macrophage', ...]

# Load species data
adata_chimp = ad.read_h5ad('chimp_data.h5ad')

# Predict with alignment
adata_pred = predict_regions(
    model=human_model,
    regions=adata_chimp,
    target_cell_types=human_celltypes,
    align_cell_types=True,
    fill_missing_with_zeros=True
)

# Compare
corr_matrix, self_corr = compare_predictions(adata_pred)

# Visualize
plt.figure(figsize=(10, 8))
sns.heatmap(corr_matrix, cmap='RdBu_r', center=0, vmin=-1, vmax=1)
plt.title(f"Chimp: Mean corr = {self_corr['correlation'].mean():.3f}")
plt.savefig('chimp_correlation.pdf')
```

---

### Workflow 3: Multi-Species Analysis

```python
from crested_utils import predict_regions, compare_predictions

species_list = ['chimp', 'gorilla', 'macaque']
results = {}

for species in species_list:
    # Load
    adata = ad.read_h5ad(f'{species}_data.h5ad')
    
    # Predict with alignment
    adata_pred = predict_regions(
        model=human_model,
        regions=adata,
        target_cell_types=human_celltypes,
        align_cell_types=True
    )
    
    # Evaluate
    corr_matrix, self_corr = compare_predictions(adata_pred)
    
    # Store
    results[species] = {
        'mean_corr': self_corr['correlation'].mean(),
        'corr_matrix': corr_matrix,
        'self_corr': self_corr
    }

# Compare across species
for species, result in results.items():
    print(f"{species}: {result['mean_corr']:.3f}")
```

---

## Tips & Tricks

### Memory Issues?

```python
# Reduce chunk size
adata_pred = predict_regions(model, adata, chunk_size=500)

# Use disk-based saving
adata_pred = predict_regions(
    model, adata, 
    chunk_size=500,
    disk_based_saving=True
)
```

### Only Want Common Cell Types?

```python
# Don't fill missing cell types
adata_pred = predict_regions(
    model, adata,
    target_cell_types=reference_celltypes,
    fill_missing_with_zeros=False  # Only keep common
)
```

### Custom Visualization?

```python
import matplotlib.pyplot as plt
import seaborn as sns

corr_matrix, self_corr = compare_predictions(adata)

# Heatmap
plt.figure(figsize=(12, 10))
sns.heatmap(
    corr_matrix,
    cmap='RdBu_r',
    center=0,
    vmin=-1,
    vmax=1,
    square=True,
    annot=True  # Add values
)

# Bar plot of self-correlations
plt.figure(figsize=(10, 6))
self_corr_sorted = self_corr.sort_values('correlation', ascending=False)
plt.barh(self_corr_sorted['cell_type'], self_corr_sorted['correlation'])
plt.xlabel('Correlation')
plt.title('Prediction Quality by Cell Type')
```

---

## Parameter Quick Guide

### chunk_size
- **1000** (default): Good for most datasets
- **2000-5000**: If you have lots of RAM
- **500**: If you get OOM errors
- **100**: For very limited memory

### batch_size
- **16** (default): Good balance
- **32-64**: Faster if you have good GPU
- **8**: If GPU memory is limited

### align_cell_types
- **True**: Align to target_cell_types (recommended for cross-species)
- **False**: Keep original cell types

### fill_missing_with_zeros
- **True**: Add missing cell types as zeros (consistent dimensions)
- **False**: Only keep common cell types (fair comparison)

### method (correlation)
- **'pearson'**: Linear relationships (default)
- **'spearman'**: Monotonic relationships (rank-based)

---

## Error Messages

### "Must provide either 'regions'"
→ Pass `regions=adata` or `regions=region_list`

### "Must provide 'genome' when regions is a list"
→ Pass `genome=crested.Genome('path/to/genome.fa')`

### "Shape mismatch"
→ Check that predictions match input dimensions

### "A and B have no columns in common"
→ DataFrames have different regions/columns

---

## Need More Help?

1. **Full docs**: See `README.md`
2. **Examples**: See `example_usage.py`
3. **Migration**: See `MIGRATION_GUIDE.md`
4. **Tests**: Run `python test_crested_utils.py`

---

## One-Liner Examples

```python
# Predict
adata_pred = predict_regions(model, regions=adata)

# Predict + align
adata_pred = predict_regions(model, regions=adata, target_cell_types=ref_cts, align_cell_types=True)

# Compare
corr, self_corr = compare_predictions(adata_pred)

# Mean correlation
mean = self_corr['correlation'].mean()

# Align separately
adata_aligned = align_adata_cell_types(adata, target_cell_types=ref_cts)
```

---

## Summary

**Main function**: `predict_regions()` - does prediction + optional alignment

**For cross-species**: Use `target_cell_types` + `align_cell_types=True`

**To evaluate**: Use `compare_predictions()` - returns matrix + statistics

**For custom work**: Use `rowwise_correlation()` directly on DataFrames

---

Keep this file handy for quick reference!
