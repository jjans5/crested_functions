# CREsted Functions# CREsted Functions# CREsted Utilities for Cross-Species Analysis



Utility functions for cross-species enhancer predictions using the [CREsted](https://github.com/aertslab/CREsted) package.



## OverviewUtility functions for cross-species enhancer predictions using the [CREsted](https://github.com/aertslab/CREsted) package.This package provides utility functions for working with [CREsted](https://github.com/aertslab/CREsted) predictions across different species, with a focus on comparing predicted vs observed values even when cell type dimensions don't match.



This repository provides tools for:

- **Predicting accessibility** of genomic regions across cell types using trained CREsted models

- **Comparing predictions** with real scATAC-seq data, even when cell types don't match## Overview## Overview

- **In-silico mutagenesis** for analyzing single-nucleotide variants

- **SNP analysis from BED files** with automated sequence extraction and testing



## InstallationThis repository provides tools for:The main module `crested_utils.py` provides:



```bash- **Predicting accessibility** of genomic regions across cell types using trained CREsted models

# Clone the repository

git clone https://github.com/jjans5/crested_functions.git- **Comparing predictions** with real scATAC-seq data, even when cell types don't match1. **`predict_regions()`** - Unified function for making predictions on regions or AnnData objects

cd crested_functions

- **In-silico mutagenesis** for analyzing single-nucleotide variants2. **`align_adata_cell_types()`** - Align cell types across species (handling missing cell types)

# Install dependencies

pip install crested anndata numpy pandas3. **`compare_predictions()`** - Compare predicted vs true values with correlation metrics

```

## Installation4. **`rowwise_correlation()`** - Compute row-wise correlations between DataFrames

## Quick Start

5. **`resize_region()`** - Resize genomic regions to specific lengths

### Option 1: Predict from Region Strings

```bash

**New!** You can now predict directly from a list of region strings:

# Clone the repository## Key Improvements Over Original Code

```python

from src.crested_utils import predict_regionsgit clone https://github.com/jjans5/crested_functions.git



# Predict from region stringscd crested_functions### 1. Unified Prediction Interface

regions = ["chr1:1000-2000", "chr1:3000-4000", "chr2:5000-6000"]

predicted_adata = predict_regions(

    model="path/to/model.keras",

    regions=regions,# Install dependencies**Before:** Multiple scattered functions with unclear interfaces

    genome=genome,

    target_cell_types=["CellTypeA", "CellTypeB"]pip install crested anndata numpy pandas

)

``````**Now:** Single `predict_regions()` function that handles:



### Option 2: Predict from AnnData- Direct AnnData objects



When your model and data have overlapping cell types:## Quick Start- Lists of region strings



```python- Automatic chunking for memory efficiency

from src.crested_utils import predict_regions, compare_predictions

### Option 1: Full Version (matching cell types)- Cell type alignment in one step

# Predict from AnnData

predicted_adata = predict_regions(

    model="path/to/model.keras",

    regions=input_adata,When your model and data have overlapping cell types:```python

    target_cell_types=["CellTypeA", "CellTypeB"]

)# Simple case: predict on existing AnnData



# Compare with real data```pythonadata_pred = predict_regions(model=model, regions=adata)

results = compare_predictions(predicted_adata, real_adata)

print(f"Mean correlation: {results['correlations'].mean():.3f}")from src.crested_utils import predict_regions, compare_predictions

```

# With cell type alignment for cross-species comparison

### Option 3: Minimal Version (zero overlap in cell types)

# Predict accessibilityadata_pred = predict_regions(

When model and data have completely different cell types:

predicted_adata = predict_regions(    model=model,

```python

from src.minimal_predict import predict_accessibility, compare_across_celltypes    model="path/to/model.keras",    regions=adata_species,



# Predict using model's own cell types    regions=input_adata,  # or just regions as strings/GenomicRanges    target_cell_types=human_celltypes,

predicted_adata = predict_accessibility(

    model="path/to/model.keras",    class_names=["CellTypeA", "CellTypeB"]  # model's cell types    align_cell_types=True

    regions=input_adata

)))



# Compare across all cell type pairs```

comparison = compare_across_celltypes(predicted_adata, real_adata)

print(f"Best match: {comparison['best_matches'].iloc[0]}")# Compare with real data

```

results = compare_predictions(predicted_adata, real_adata)### 2. Automatic Cell Type Alignment

### In-Silico Mutagenesis

print(f"Mean correlation: {results['correlations'].mean():.3f}")

**Updated with log2 fold changes!**

```**Before:** Manual, error-prone code with many steps to align cell types

Analyze effects of single-nucleotide variants:



```python

from src.insilico_mutagenesis_vect import insilico_mutagenesis_vect### Option 2: Minimal Version (zero overlap in cell types)**Now:** Automatic alignment with options:



results = insilico_mutagenesis_vect(- Fill missing cell types with zeros (for comparison)

    seq="ACGT...",  # Your sequence

    model=model,When model and data have completely different cell types:- Or keep only common cell types

    adata=adata,

    chrom="chr1",- Preserves all layers (predictions, weights, etc.)

    start=1000

)```python- Maintains exact order specified



# Returns wildtype predictions and all possible mutationsfrom src.minimal_predict import predict_accessibility, compare_across_celltypes

print(results['wildtype'])  # Original predictions

print(results['mut_df_long'])  # Includes 'log2fc' and 'delta' columns```python

```

# Predict using model's own cell types# Align species data to human cell type order

### SNP Analysis from BED Files

predicted_adata = predict_accessibility(adata_aligned = align_adata_cell_types(

**New function!** Analyze SNPs directly from BED files:

    model="path/to/model.keras",    adata_species,

```python

from src.insilico_mutagenesis_vect import snp_mutagenesis_from_bed    regions=input_adata    target_cell_types=human_celltypes,



# Analyze all SNPs in a BED file)    fill_missing=True  # Add zeros for missing cell types

results = snp_mutagenesis_from_bed(

    bed_file="snps.bed",  # Can have just positions or include ref/alt)

    model=model,

    adata=adata,# Compare across all cell type pairs```

    genome=genome,

    seq_length=2114  # Optional, inferred from model if not providedcomparison = compare_across_celltypes(predicted_adata, real_adata)

)

print(f"Best match: {comparison['best_matches'].iloc[0]}")### 3. Built-in Comparison Functions

# Results include log2fc for each SNP-cell_type combination

print(results.head())```

```

**Before:** Manual DataFrame creation and correlation calculation

**BED file formats supported:**

- **3 columns:** `chrom start end` - will test all possible mutations at each position### In-Silico Mutagenesis

- **5 columns:** `chrom start end ref alt` - will test specific ref>alt changes

- **6+ columns:** `chrom start end name ref alt` - additional columns are ignored**Now:** Direct comparison with built-in functions:



## Repository StructureAnalyze effects of single-nucleotide variants:



``````python

crested_functions/

├── src/                           # Source code```python# Compare predictions vs true values

│   ├── crested_utils.py          # Full version with alignment

│   ├── minimal_predict.py        # Minimal version for zero overlapfrom src.insilico_mutagenesis_vect import insilico_mutagenesis_vectcorr_matrix, self_corr = compare_predictions(

│   └── insilico_mutagenesis_vect.py  # SNP analysis with log2fc

├── tests/                         # Unit tests    adata,

│   └── test_crested_utils.py     

├── scripts/                       # Example scriptsresults = insilico_mutagenesis_vect(    prediction_layer="predicted",

│   ├── example_usage.py          # Full version examples

│   ├── demo_minimal.py           # Minimal version demo    model="path/to/model.keras",    method="pearson"

│   └── insilico_mutagenesis_vect_example.py

├── LICENSE                        # MIT License    region="chr1:1000-2000",)

└── README.md                      # This file

```    cell_types=["CellTypeA", "CellTypeB"]



## Key Features)print(f"Mean self-correlation: {self_corr['correlation'].mean():.3f}")



### Enhanced predict_regions()```

- ✅ **New:** Accept list of region strings directly

- ✅ Accept AnnData objects# Returns wildtype predictions and all possible mutations

- ✅ Automatic memory-efficient chunking

- ✅ Cell type alignment in one stepprint(results['wildtype'])  # Original predictions### 4. Improved Memory Management

- ✅ Optional genome parameter for region strings

print(results['mutations_wide'])  # All mutations in wide format

### Improved In-Silico Mutagenesis

- ✅ **New:** Log2 fold change calculation (log2(mutant/wildtype))```- Automatic memory monitoring (if `psutil` available)

- ✅ Absolute delta values for reference

- ✅ Vectorized predictions for speed- Adaptive chunk sizing based on available memory

- ✅ Genomic position tracking

## Repository Structure- Option for disk-based processing for very large datasets

### New SNP Analysis Function

- ✅ **New:** Process BED files with SNP positions- Progress bars for long operations

- ✅ Automatic sequence extraction centered on SNPs

- ✅ Support for ref/alt or test all mutations```

- ✅ Returns log2fc for each cell type

- ✅ Flexible BED format parsingcrested_functions/### 5. Better Error Handling and Logging



## Key Functions├── src/                           # Source code



### crested_utils.py│   ├── crested_utils.py          # Full version with alignment- Comprehensive logging at each step

- `predict_regions()`: Predict accessibility for genomic regions (now accepts region lists!)

- `align_adata_cell_types()`: Align cell types between datasets│   ├── minimal_predict.py        # Minimal version for zero overlap- Clear error messages

- `compare_predictions()`: Compare predicted vs real values

- `rowwise_correlation()`: Compute row-wise correlations│   └── insilico_mutagenesis_vect.py  # SNP analysis- Shape validation



### minimal_predict.py├── tests/                         # Unit tests- Informative progress bars

- `predict_accessibility()`: Predict using model's cell types

- `compare_across_celltypes()`: All-vs-all cell type comparison│   └── test_crested_utils.py     

- `find_best_matches()`: Identify best matching cell types

├── scripts/                       # Example scripts### 6. Type Hints and Documentation

### insilico_mutagenesis_vect.py

- `insilico_mutagenesis_vect()`: Vectorized SNP effect analysis with log2fc│   ├── example_usage.py          # Full version examples

- `snp_mutagenesis_from_bed()`: Process SNPs from BED files (NEW!)

│   ├── demo_minimal.py           # Minimal version demo- Full type hints for better IDE support

## Which Version Should I Use?

│   └── insilico_mutagenesis_vect_example.py- Comprehensive docstrings with examples

- **Use `crested_utils.py`** if your model and data share some cell types (e.g., both have "CD4 T cells")

- **Use `minimal_predict.py`** if your model was trained on one species/dataset and you're predicting on a completely different species/dataset with no cell type overlap├── LICENSE                        # MIT License- Clear parameter descriptions

- **Use `snp_mutagenesis_from_bed()`** if you have a BED file with SNP positions and want to analyze their effects

└── README.md                      # This file

## What's New

```## Installation

### Version 1.1.0 (Latest)

- ✨ `predict_regions()` now accepts list of region strings as input

- ✨ New function: `snp_mutagenesis_from_bed()` for BED file SNP analysis

- ✨ In-silico mutagenesis now returns log2 fold changes (log2fc)## Key Functions```bash

- ✨ Improved documentation and examples

# Install CREsted (follow their instructions)

## Requirements

### crested_utils.pypip install crested

- Python 3.9+

- crested- `predict_regions()`: Predict accessibility for genomic regions

- anndata

- numpy- `align_adata_cell_types()`: Align cell types between datasets# Optional: for memory monitoring

- pandas

- psutil (optional, for memory monitoring)- `compare_predictions()`: Compare predicted vs real valuespip install psutil



## License- `rowwise_correlation()`: Compute row-wise correlations```



MIT License - see [LICENSE](LICENSE) file for details.



## Citation### minimal_predict.py## Quick Start



If you use these functions, please cite the CREsted paper:- `predict_accessibility()`: Predict using model's cell types

- CREsted: [https://github.com/aertslab/CREsted](https://github.com/aertslab/CREsted)

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
