# CREsted Functions# CREsted Utilities for Cross-Species Analysis



Utility functions for cross-species enhancer predictions using the [CREsted](https://github.com/aertslab/CREsted) package.This package provides utility functions for working with [CREsted](https://github.com/aertslab/CREsted) predictions across different species, with a focus on comparing predicted vs observed values even when cell type dimensions don't match.



## Overview## Overview



This repository provides tools for:The main module `crested_utils.py` provides:

- **Predicting accessibility** of genomic regions across cell types using trained CREsted models

- **Comparing predictions** with real scATAC-seq data, even when cell types don't match1. **`predict_regions()`** - Unified function for making predictions on regions or AnnData objects

- **In-silico mutagenesis** for analyzing single-nucleotide variants2. **`align_adata_cell_types()`** - Align cell types across species (handling missing cell types)

3. **`compare_predictions()`** - Compare predicted vs true values with correlation metrics

## Installation4. **`rowwise_correlation()`** - Compute row-wise correlations between DataFrames

5. **`resize_region()`** - Resize genomic regions to specific lengths

```bash

# Clone the repository## Key Improvements Over Original Code

git clone https://github.com/jjans5/crested_functions.git

cd crested_functions### 1. Unified Prediction Interface



# Install dependencies**Before:** Multiple scattered functions with unclear interfaces

pip install crested anndata numpy pandas

```**Now:** Single `predict_regions()` function that handles:

- Direct AnnData objects

## Quick Start- Lists of region strings

- Automatic chunking for memory efficiency

### Option 1: Full Version (matching cell types)- Cell type alignment in one step



When your model and data have overlapping cell types:```python

# Simple case: predict on existing AnnData

```pythonadata_pred = predict_regions(model=model, regions=adata)

from src.crested_utils import predict_regions, compare_predictions

# With cell type alignment for cross-species comparison

# Predict accessibilityadata_pred = predict_regions(

predicted_adata = predict_regions(    model=model,

    model="path/to/model.keras",    regions=adata_species,

    regions=input_adata,  # or just regions as strings/GenomicRanges    target_cell_types=human_celltypes,

    class_names=["CellTypeA", "CellTypeB"]  # model's cell types    align_cell_types=True

))

```

# Compare with real data

results = compare_predictions(predicted_adata, real_adata)### 2. Automatic Cell Type Alignment

print(f"Mean correlation: {results['correlations'].mean():.3f}")

```**Before:** Manual, error-prone code with many steps to align cell types



### Option 2: Minimal Version (zero overlap in cell types)**Now:** Automatic alignment with options:

- Fill missing cell types with zeros (for comparison)

When model and data have completely different cell types:- Or keep only common cell types

- Preserves all layers (predictions, weights, etc.)

```python- Maintains exact order specified

from src.minimal_predict import predict_accessibility, compare_across_celltypes

```python

# Predict using model's own cell types# Align species data to human cell type order

predicted_adata = predict_accessibility(adata_aligned = align_adata_cell_types(

    model="path/to/model.keras",    adata_species,

    regions=input_adata    target_cell_types=human_celltypes,

)    fill_missing=True  # Add zeros for missing cell types

)

# Compare across all cell type pairs```

comparison = compare_across_celltypes(predicted_adata, real_adata)

print(f"Best match: {comparison['best_matches'].iloc[0]}")### 3. Built-in Comparison Functions

```

**Before:** Manual DataFrame creation and correlation calculation

### In-Silico Mutagenesis

**Now:** Direct comparison with built-in functions:

Analyze effects of single-nucleotide variants:

```python

```python# Compare predictions vs true values

from src.insilico_mutagenesis_vect import insilico_mutagenesis_vectcorr_matrix, self_corr = compare_predictions(

    adata,

results = insilico_mutagenesis_vect(    prediction_layer="predicted",

    model="path/to/model.keras",    method="pearson"

    region="chr1:1000-2000",)

    cell_types=["CellTypeA", "CellTypeB"]

)print(f"Mean self-correlation: {self_corr['correlation'].mean():.3f}")

```

# Returns wildtype predictions and all possible mutations

print(results['wildtype'])  # Original predictions### 4. Improved Memory Management

print(results['mutations_wide'])  # All mutations in wide format

```- Automatic memory monitoring (if `psutil` available)

- Adaptive chunk sizing based on available memory

## Repository Structure- Option for disk-based processing for very large datasets

- Progress bars for long operations

```

crested_functions/### 5. Better Error Handling and Logging

├── src/                           # Source code

│   ├── crested_utils.py          # Full version with alignment- Comprehensive logging at each step

│   ├── minimal_predict.py        # Minimal version for zero overlap- Clear error messages

│   └── insilico_mutagenesis_vect.py  # SNP analysis- Shape validation

├── tests/                         # Unit tests- Informative progress bars

│   └── test_crested_utils.py     

├── scripts/                       # Example scripts### 6. Type Hints and Documentation

│   ├── example_usage.py          # Full version examples

│   ├── demo_minimal.py           # Minimal version demo- Full type hints for better IDE support

│   └── insilico_mutagenesis_vect_example.py- Comprehensive docstrings with examples

├── LICENSE                        # MIT License- Clear parameter descriptions

└── README.md                      # This file

```## Installation



## Key Functions```bash

# Install CREsted (follow their instructions)

### crested_utils.pypip install crested

- `predict_regions()`: Predict accessibility for genomic regions

- `align_adata_cell_types()`: Align cell types between datasets# Optional: for memory monitoring

- `compare_predictions()`: Compare predicted vs real valuespip install psutil

- `rowwise_correlation()`: Compute row-wise correlations```



### minimal_predict.py## Quick Start

- `predict_accessibility()`: Predict using model's cell types

- `compare_across_celltypes()`: All-vs-all cell type comparison### Basic Prediction

- `find_best_matches()`: Identify best matching cell types

```python

### insilico_mutagenesis_vect.pyfrom crested_utils import predict_regions

- `insilico_mutagenesis_vect()`: Vectorized SNP effect analysisimport crested

import anndata as ad

## Which Version Should I Use?

# Load data and model

- **Use `crested_utils.py`** if your model and data share some cell types (e.g., both have "CD4 T cells")adata = ad.read_h5ad("your_data.h5ad")

- **Use `minimal_predict.py`** if your model was trained on one species/dataset and you're predicting on a completely different species/dataset with no cell type overlapmodel = crested.load_model("your_model")



## Requirements# Make predictions (with automatic chunking)

adata_pred = predict_regions(

- Python 3.9+    model=model,

- crested    regions=adata,

- anndata    chunk_size=1000,  # Adjust based on your memory

- numpy    verbose=True

- pandas)

- psutil (optional, for memory monitoring)

# Access predictions

## Licenseprint(adata_pred.layers['predicted'].shape)

```

MIT License - see [LICENSE](LICENSE) file for details.

### Cross-Species Comparison

## Citation

```python

If you use these functions, please cite the CREsted paper:from crested_utils import predict_regions, compare_predictions

- CREsted: [https://github.com/aertslab/CREsted](https://github.com/aertslab/CREsted)import matplotlib.pyplot as plt

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
