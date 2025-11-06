# CREsted Utilities for Cross-Species Analysis

This package provides utility functions for working with [CREsted](https://github.com/aertslab/CREsted) predictions across different species, with a focus on comparing predicted vs observed values even when cell type dimensions don't match.

## Overview

The main module `crested_utils.py` provides:

1. **`predict_regions()`** - Unified function for making predictions on regions or AnnData objects
2. **`align_adata_cell_types()`** - Align cell types across species (handling missing cell types)
3. **`compare_predictions()`** - Compare predicted vs true values with correlation metrics
4. **`rowwise_correlation()`** - Compute row-wise correlations between DataFrames
5. **`resize_region()`** - Resize genomic regions to specific lengths

## Key Improvements Over Original Code

### 1. Unified Prediction Interface

**Before:** Multiple scattered functions with unclear interfaces

**Now:** Single `predict_regions()` function that handles:
- Direct AnnData objects
- Lists of region strings
- Automatic chunking for memory efficiency
- Cell type alignment in one step

```python
# Simple case: predict on existing AnnData
adata_pred = predict_regions(model=model, regions=adata)

# With cell type alignment for cross-species comparison
adata_pred = predict_regions(
    model=model,
    regions=adata_species,
    target_cell_types=human_celltypes,
    align_cell_types=True
)
```

### 2. Automatic Cell Type Alignment

**Before:** Manual, error-prone code with many steps to align cell types

**Now:** Automatic alignment with options:
- Fill missing cell types with zeros (for comparison)
- Or keep only common cell types
- Preserves all layers (predictions, weights, etc.)
- Maintains exact order specified

```python
# Align species data to human cell type order
adata_aligned = align_adata_cell_types(
    adata_species,
    target_cell_types=human_celltypes,
    fill_missing=True  # Add zeros for missing cell types
)
```

### 3. Built-in Comparison Functions

**Before:** Manual DataFrame creation and correlation calculation

**Now:** Direct comparison with built-in functions:

```python
# Compare predictions vs true values
corr_matrix, self_corr = compare_predictions(
    adata,
    prediction_layer="predicted",
    method="pearson"
)

print(f"Mean self-correlation: {self_corr['correlation'].mean():.3f}")
```

### 4. Improved Memory Management

- Automatic memory monitoring (if `psutil` available)
- Adaptive chunk sizing based on available memory
- Option for disk-based processing for very large datasets
- Progress bars for long operations

### 5. Better Error Handling and Logging

- Comprehensive logging at each step
- Clear error messages
- Shape validation
- Informative progress bars

### 6. Type Hints and Documentation

- Full type hints for better IDE support
- Comprehensive docstrings with examples
- Clear parameter descriptions

## Installation

```bash
# Install CREsted (follow their instructions)
pip install crested

# Optional: for memory monitoring
pip install psutil
```

## Quick Start

### Basic Prediction

```python
from crested_utils import predict_regions
import crested
import anndata as ad

# Load data and model
adata = ad.read_h5ad("your_data.h5ad")
model = crested.load_model("your_model")

# Make predictions (with automatic chunking)
adata_pred = predict_regions(
    model=model,
    regions=adata,
    chunk_size=1000,  # Adjust based on your memory
    verbose=True
)

# Access predictions
print(adata_pred.layers['predicted'].shape)
```

### Cross-Species Comparison

```python
from crested_utils import predict_regions, compare_predictions
import matplotlib.pyplot as plt
import seaborn as sns

# Define reference cell types (e.g., from human)
human_celltypes = ['T_cell', 'B_cell', 'NK_cell', 'Macrophage']

# Load species data
adata_chimp = ad.read_h5ad("chimpanzee_data.h5ad")

# Predict with cell type alignment
adata_pred = predict_regions(
    model=model,
    regions=adata_chimp,
    target_cell_types=human_celltypes,
    align_cell_types=True,
    fill_missing_with_zeros=True
)

# Compare predictions
corr_matrix, self_corr = compare_predictions(adata_pred)

# Visualize
plt.figure(figsize=(10, 8))
sns.heatmap(corr_matrix, cmap='RdBu_r', center=0, vmin=-1, vmax=1)
plt.title('Predicted vs True Correlation')
plt.savefig('correlation.pdf')
```

## Complete Examples

See `example_usage.py` for comprehensive examples including:

1. Basic prediction on AnnData
2. Cross-species prediction with alignment
3. Comparing predictions vs true values
4. Full cross-species analysis pipeline
5. Manual cell type alignment
6. Region resizing
7. Direct row-wise correlations

## Function Reference

### `predict_regions()`

Main function for making predictions on genomic regions.

**Parameters:**
- `model`: CREsted model
- `regions`: List of region strings OR AnnData object
- `genome`: crested.Genome (required if regions is a list)
- `batch_size`: Batch size for prediction (default: 16)
- `chunk_size`: Regions per chunk (default: 1000)
- `layer_name`: Output layer name (default: 'predicted')
- `target_cell_types`: Target cell types for alignment
- `align_cell_types`: Whether to align (default: True)
- `fill_missing_with_zeros`: Fill missing cell types (default: True)
- `disk_based_saving`: Use disk for large datasets (default: False)
- `verbose`: Enable logging (default: True)

**Returns:** AnnData with predictions in `.layers[layer_name]`

### `align_adata_cell_types()`

Align AnnData cell types to a target list.

**Parameters:**
- `adata`: Input AnnData
- `target_cell_types`: Desired cell types in order
- `fill_missing`: Add zeros for missing (default: True)
- `verbose`: Enable logging (default: True)

**Returns:** New AnnData with aligned cell types

### `compare_predictions()`

Compare predicted vs true values.

**Parameters:**
- `adata`: AnnData with predictions
- `prediction_layer`: Name of prediction layer (default: 'predicted')
- `true_layer`: Name of true layer or None for .X (default: None)
- `cell_types`: Cell types to compare or None for all (default: None)
- `method`: 'pearson' or 'spearman' (default: 'pearson')

**Returns:** 
- `corr_df`: Full correlation matrix
- `diag_corr`: DataFrame with diagonal (self) correlations

### `rowwise_correlation()`

Compute row-wise correlation between two DataFrames.

**Parameters:**
- `A`: First DataFrame (n_rows_A x n_columns)
- `B`: Second DataFrame (n_rows_B x n_columns)
- `method`: 'pearson' or 'spearman' (default: 'pearson')

**Returns:** Correlation matrix (n_rows_A x n_rows_B)

### `resize_region()`

Resize a genomic region string.

**Parameters:**
- `region`: Region string (e.g., 'chr1:1000-2000')
- `new_length`: Desired length
- `summit`: Center position or None for midpoint (default: None)

**Returns:** Resized region string

## Use Cases

### 1. Cross-Species Prediction

When comparing regulatory predictions across species where not all cell types are present:

```python
# Each species may have different cell types
# This ensures all predictions have the same dimensions
for species in ['chimp', 'gorilla', 'macaque']:
    adata_species = ad.read_h5ad(f"{species}_data.h5ad")
    adata_pred = predict_regions(
        model=human_model,
        regions=adata_species,
        target_cell_types=human_celltypes,
        align_cell_types=True
    )
    # Now all have same shape for comparison
```

### 2. Model Comparison

Compare predictions from different models:

```python
# Get predictions from two models
adata_pred1 = predict_regions(model1, regions=adata)
adata_pred2 = predict_regions(model2, regions=adata)

# Compare them
pred1_df = pd.DataFrame(adata_pred1.layers['predicted'], index=adata.obs_names)
pred2_df = pd.DataFrame(adata_pred2.layers['predicted'], index=adata.obs_names)
corr = rowwise_correlation(pred1_df, pred2_df)
```

### 3. Large Dataset Processing

For very large datasets that don't fit in memory:

```python
adata_pred = predict_regions(
    model=model,
    regions=large_adata,
    chunk_size=500,  # Smaller chunks
    disk_based_saving=True,  # Save chunks to disk
    temp_dir="/path/to/temp"  # Specify temp location
)
```

## Comparison with Original Code

### Original Code Issues

1. **Scattered functions**: `predict_chunked.py`, `species_pred.py`, `rowwise_corr.py` were separate
2. **Manual alignment**: Required ~30 lines of code to align cell types
3. **No validation**: Shape mismatches caused silent errors
4. **Poor documentation**: Unclear function purposes and parameters
5. **Error-prone**: Easy to make mistakes with manual DataFrame manipulation

### New Code Benefits

1. **Unified interface**: Single import, clear API
2. **Automatic alignment**: One parameter to align cell types
3. **Built-in validation**: Checks shapes, provides clear errors
4. **Comprehensive docs**: Every function documented with examples
5. **Type-safe**: Type hints catch errors early

## Tips and Best Practices

1. **Memory**: Start with `chunk_size=1000`, reduce if you get OOM errors
2. **Cell types**: Always specify `target_cell_types` for cross-species work
3. **Validation**: Check correlations - low values indicate issues
4. **Logging**: Keep `verbose=True` initially to monitor progress
5. **Disk-based**: Use for datasets > 50k regions or low memory systems

## Backward Compatibility

The old function names are aliased for backward compatibility:

```python
# Old names still work
from crested_utils import predict_chunked, rowwise_corr

# But use new names in new code
from crested_utils import predict_regions, rowwise_correlation
```

## License

Same as CREsted package.

## Citation

If you use these utilities, please cite the CREsted package:

> Kempynck, N., De Winter, S., et al. CREsted: modeling genomic and synthetic cell type-specific enhancers across tissues and species. bioRxiv (2025).

## Support

For issues with these utilities, open an issue in your repository.
For issues with CREsted itself, see: https://github.com/aertslab/CREsted/issues
