# Which Version to Use?

## Decision Tree

```
Do your model and data share cell types?
│
├─ YES → Use crested_utils.py
│         • predict_regions() with align_cell_types=True
│         • Automatic alignment to common cell types
│         • Example: Human model on chimp data
│
└─ NO  → Use minimal_predict.py
          • predict_accessibility()
          • compare_across_celltypes()
          • Example: Mouse model on zebrafish data
```

## crested_utils.py

**Use when:**
- Model and data have some overlapping cell types
- Cross-species comparisons (human ↔ chimp)
- You want automatic cell type alignment

**Example:**
```python
from src.crested_utils import predict_regions

adata_pred = predict_regions(
    model=human_model,
    regions=chimp_data,
    target_cell_types=["T_cell", "B_cell", "NK_cell"],
    align_cell_types=True
)
```

## minimal_predict.py

**Use when:**
- Model and data have NO overlapping cell types
- Different species with different nomenclature
- You want to find functional equivalents

**Example:**
```python
from src.minimal_predict import predict_accessibility, compare_across_celltypes

# Predict with model's cell types
adata_pred = predict_accessibility(model, data)

# Find best matches across all pairs
comparison = compare_across_celltypes(adata_pred, adata_real)
```

## Key Differences

| Feature | crested_utils | minimal_predict |
|---------|--------------|-----------------|
| Cell type alignment | Automatic | Manual/all-pairs |
| Input format | Regions or AnnData | AnnData only |
| Output cell types | Aligned to target | Model's native |
| Use case | Overlapping types | Zero overlap |
| Complexity | Higher | Simpler |
