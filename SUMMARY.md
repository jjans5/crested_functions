# Summary of Improvements to CREsted Utilities

## What Was Done

I reviewed your scattered CREsted utility functions and consolidated them into a unified, well-documented module with significant improvements.

## Files Created

1. **`crested_utils.py`** - Main utility module with all functions
2. **`example_usage.py`** - 7 comprehensive examples
3. **`README.md`** - Complete documentation
4. **`MIGRATION_GUIDE.md`** - Step-by-step migration guide
5. **`SUMMARY.md`** - This file

## Key Improvements

### 1. Code Consolidation
- **Before**: 4 separate files with unclear relationships
- **After**: 1 unified module with clear API

### 2. Reduced Code Complexity
- **Cell type alignment**: 30+ lines → 4 lines (87% reduction)
- **Complete pipeline**: ~100 lines → ~40 lines (60% reduction)

### 3. New Main Function: `predict_regions()`

This is your **primary interface** for predictions:

```python
from crested_utils import predict_regions

# Simple prediction
adata_pred = predict_regions(model, regions=adata)

# With cell type alignment for cross-species
adata_pred = predict_regions(
    model=model,
    regions=adata_species,
    target_cell_types=human_celltypes,
    align_cell_types=True,
    fill_missing_with_zeros=True
)
```

**Handles**:
- AnnData objects or region lists
- Automatic chunking for large datasets
- Cell type alignment in one step
- Memory monitoring and optimization
- Progress tracking
- Disk-based processing option

### 4. Cell Type Alignment Made Easy

**Your original manual alignment code** (from `species_pred.py`):
- 30+ lines of error-prone NumPy/Pandas manipulation
- Manual handling of missing cell types
- Manual preservation of layers and obsm
- Easy to make mistakes

**New function**:
```python
adata_aligned = align_adata_cell_types(
    adata_species,
    target_cell_types=human_celltypes,
    fill_missing=True
)
```

**Automatically**:
- ✅ Identifies present and missing cell types
- ✅ Fills missing with zeros (or excludes them)
- ✅ Reorders to exact target order
- ✅ Preserves all layers (predictions, etc.)
- ✅ Preserves obsm (weights, etc.)
- ✅ Validates dtypes and shapes

### 5. Built-in Comparison Function

**Before**: Manual DataFrame creation + separate correlation function

**Now**: Single function with statistics

```python
corr_matrix, self_corr = compare_predictions(
    adata,
    prediction_layer='predicted',
    cell_types=common_celltypes
)

print(f"Mean correlation: {self_corr['correlation'].mean():.3f}")
print(self_corr.sort_values('correlation'))
```

Returns:
- Full correlation matrix (for heatmaps)
- Diagonal correlations (how well each cell type is predicted)
- Already formatted for visualization

### 6. Better Memory Management

Your original `predict_chunked()` had good memory features. The new version adds:
- Automatic chunk size adjustment based on available memory
- Better progress reporting
- More efficient concatenation
- Clearer logging at each step

### 7. Improved Developer Experience

**Type Hints**:
```python
def predict_regions(
    model,
    regions: Optional[Union[List[str], ad.AnnData]] = None,
    target_cell_types: Optional[List[str]] = None,
    align_cell_types: bool = True,
    verbose: bool = True,
) -> ad.AnnData:
```

**Comprehensive Docstrings**:
- Parameter descriptions
- Return value descriptions
- Usage examples
- Cross-references

**Better Error Messages**:
```python
if regions is None:
    raise ValueError("Must provide either 'regions' (list or AnnData)")

if adata_concat.shape != adata.shape:
    raise ValueError(f"Shape mismatch: got {adata_concat.shape} expected {adata.shape}")
```

## Function Summary

### Core Functions

| Function | Purpose | Key Feature |
|----------|---------|-------------|
| `predict_regions()` | Make predictions on regions/AnnData | Unified interface, auto-alignment |
| `align_adata_cell_types()` | Align cell types across datasets | Handles missing cell types |
| `compare_predictions()` | Compare predicted vs true | Returns matrix + statistics |
| `rowwise_correlation()` | Row-wise correlation | Fast vectorized implementation |
| `resize_region()` | Resize genomic regions | Unchanged from original |

### Internal Functions

- `_predict_chunked()`: Memory-efficient chunked prediction
  (Your original with improvements, now internal)

## Usage Patterns

### Pattern 1: Basic Prediction
```python
adata_pred = predict_regions(model, regions=adata)
```

### Pattern 2: Cross-Species (YOUR MAIN USE CASE)
```python
# This is what you want!
adata_pred = predict_regions(
    model=human_model,
    regions=adata_species,
    target_cell_types=human_celltypes,  # Reference cell types
    align_cell_types=True,              # Align to reference
    fill_missing_with_zeros=True        # Fill missing with zeros
)

# Now all species have same dimensions
# Easy to compare predictions vs true values
corr_matrix, self_corr = compare_predictions(adata_pred)
```

### Pattern 3: Evaluation
```python
corr_matrix, self_corr = compare_predictions(adata)
mean_corr = self_corr['correlation'].mean()
print(f"Model performance: {mean_corr:.3f}")
```

## Solving Your Original Goal

> "The end goal is to compare predicted with real values, even if the dimensions of number of cell types don't match"

**Solution**: Use `predict_regions()` with alignment

```python
# Define reference cell types (e.g., from human)
reference_celltypes = ['T_cell', 'B_cell', 'NK_cell', 'Macrophage', ...]

# For each species
for species in ['chimpanzee', 'gorilla', 'macaque']:
    adata_species = load_species_data(species)
    
    # Predict with alignment - handles missing cell types automatically
    adata_pred = predict_regions(
        model=human_model,
        regions=adata_species,
        target_cell_types=reference_celltypes,
        align_cell_types=True,
        fill_missing_with_zeros=True  # Species-specific cell types → zeros
    )
    
    # All species now have same shape: (n_reference_celltypes, n_regions)
    # Missing cell types are filled with zeros
    # Easy to compare!
    
    corr_matrix, self_corr = compare_predictions(
        adata_pred,
        cell_types=[ct for ct in reference_celltypes if ct in adata_species.obs_names]
    )
    
    print(f"{species}: mean corr = {self_corr['correlation'].mean():.3f}")
```

## Backward Compatibility

Old function names still work:
```python
# Old names (still work)
from crested_utils import predict_chunked, rowwise_corr

# New names (recommended)
from crested_utils import predict_regions, rowwise_correlation
```

## Documentation Provided

1. **README.md**: Full documentation with examples
2. **example_usage.py**: 7 complete examples covering all use cases
3. **MIGRATION_GUIDE.md**: Step-by-step migration from old to new code
4. **Inline docstrings**: Every function fully documented

## Next Steps

### Immediate Use

Start using right away:

```python
from crested_utils import predict_regions, compare_predictions

# Your main workflow
adata_pred = predict_regions(
    model, 
    regions=adata_species,
    target_cell_types=human_celltypes,
    align_cell_types=True
)

corr_matrix, self_corr = compare_predictions(adata_pred)
```

### Testing

1. Run on small dataset first
2. Compare with your old code results
3. Verify shapes and correlations match
4. Then use for full analysis

### Migration

1. Start with one species/dataset
2. Use the MIGRATION_GUIDE.md
3. Keep old files until verified
4. Migrate gradually

## Benefits Summary

✅ **60% less code** for typical workflows
✅ **Automatic cell type alignment** (no manual code needed)
✅ **Built-in validation** (catches errors early)
✅ **Better error messages** (know what went wrong)
✅ **Type hints** (IDE autocomplete & type checking)
✅ **Progress bars** (know how long it takes)
✅ **Memory optimization** (auto-adjust chunk sizes)
✅ **Comprehensive docs** (examples for everything)
✅ **Unified interface** (one import, clear API)
✅ **Backward compatible** (old code still works)

## Questions to Consider

1. **Should I use `fill_missing_with_zeros=True` or `False`?**
   - `True`: Keep all reference cell types, fill missing with zeros (good for consistent dimensions)
   - `False`: Only keep common cell types (good for fair comparison)

2. **What chunk_size should I use?**
   - Start with 1000
   - Reduce if you get out-of-memory errors
   - Increase if you have lots of RAM and want faster processing

3. **Should I use `disk_based_saving`?**
   - No: For datasets < 50k regions
   - Yes: For very large datasets or limited RAM

4. **What correlation method?**
   - `pearson`: For linear relationships (default)
   - `spearman`: For monotonic but non-linear relationships

## File Organization

Your workspace now has:

```
crested_mod/
├── crested_utils.py          # Main module (USE THIS)
├── example_usage.py          # 7 complete examples
├── README.md                 # Full documentation
├── MIGRATION_GUIDE.md        # Migration guide
├── SUMMARY.md                # This file
├── predict_chunked.py        # Original (keep for reference)
├── species_pred.py           # Original (keep for reference)
├── rowwise_corr.py          # Original (keep for reference)
└── resize_region.py         # Original (keep for reference)
```

**Recommendation**: Use `crested_utils.py` for all new code. Keep old files for reference until you've verified everything works.

---

## Quick Start

```python
# Import
from crested_utils import predict_regions, compare_predictions
import crested
import anndata as ad

# Load
model = crested.load_model("your_model")
adata = ad.read_h5ad("your_data.h5ad")

# Predict (with alignment if needed)
adata_pred = predict_regions(
    model=model,
    regions=adata,
    target_cell_types=reference_celltypes,  # Optional
    align_cell_types=True,                   # Optional
)

# Evaluate
corr_matrix, self_corr = compare_predictions(adata_pred)
print(f"Performance: {self_corr['correlation'].mean():.3f}")
```

That's it! Your original scattered functions are now a professional, well-documented utility module.
