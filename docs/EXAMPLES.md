# Examples

## Basic Prediction

```python
from src.crested_utils import predict_regions
import crested

# Load model and data
model = crested.load_model("model.keras")
adata = ad.read_h5ad("data.h5ad")

# Make predictions
adata_pred = predict_regions(model, adata)
```

## Predict from Region Strings

```python
from src.crested_utils import predict_regions

regions = ["chr1:1000-2000", "chr1:3000-4000"]
adata = predict_regions(
    model=model,
    regions=regions,
    genome=genome,
    target_cell_types=["T_cell", "B_cell"]
)
```

## Cross-Species Comparison

```python
from src.crested_utils import predict_regions, compare_predictions

# Predict with cell type alignment
adata_pred = predict_regions(
    model=model,
    regions=adata_species,
    target_cell_types=["T_cell", "B_cell", "NK_cell"],
    align_cell_types=True
)

# Compare with real data
corr_matrix, self_corr = compare_predictions(adata_pred, adata_real)
print(f"Mean correlation: {self_corr['correlation'].mean():.3f}")
```

## Zero Cell Type Overlap

When model and data have completely different cell types:

```python
from src.minimal_predict import predict_accessibility, compare_across_celltypes

# Predict using model's own cell types
adata_pred = predict_accessibility(model, adata)

# Compare across all cell type pairs
comparison = compare_across_celltypes(adata_pred, adata_real)
best_matches = comparison['best_matches']
```

## In-Silico Mutagenesis

```python
from src.insilico_mutagenesis_vect import insilico_mutagenesis_vect

# Single sequence analysis
results = insilico_mutagenesis_vect(
    seq="ACGTACGT...",
    model=model,
    adata=adata,
    chrom="chr1",
    start=1000
)

# Check results
print(results['wt_pred'])          # Wildtype predictions
print(results['mut_df_wide'])      # All mutations (wide format)
print(results['mut_df_long'])      # Tidy format with log2fc and delta
```

## SNP Analysis from BED File

```python
from src.insilico_mutagenesis_vect import snp_mutagenesis_from_bed

# Analyze SNPs
results = snp_mutagenesis_from_bed(
    bed_file="snps.bed",
    model=model,
    adata=adata,
    genome=genome
)

# Find significant effects
significant = results[abs(results['log2fc']) > 1]
print(f"Found {len(significant)} significant mutations")

# Top effects per cell type
for ct in results['cell_type'].unique():
    ct_results = results[results['cell_type'] == ct]
    top = ct_results.nlargest(5, 'log2fc')
    print(f"\n{ct} top effects:")
    print(top[['mut_id', 'log2fc']])
```

## BED File Formats

### Format 1: Positions only (tests all mutations)
```
chr1    1000    1001
chr1    2000    2001
```

### Format 2: With ref/alt (tests specific mutation)
```
chr1    1000    1001    A    G
chr1    2000    2001    C    T
```

### Format 3: With names
```
chr1    1000    1001    rs123    A    G
chr2    3000    3001    rs456    C    T
```
