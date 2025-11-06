# When to Use Which Version

Quick guide to choose between minimal and full versions.

## Decision Tree

```
Do your cell types match between model and data?
│
├─ NO (zero overlap, e.g., cross-species) ───► Use MINIMAL version
│   └─ Examples:
│      - Human model, chimp data
│      - Different cell type annotations
│      - Want to find functional matches
│
└─ YES (or can be aligned) ───────────────────► Use FULL version
    └─ Examples:
       - Same species, same annotations
       - Can map cell types (e.g., both have "T_cell")
       - Want predictions in same dimensions as input
```

---

## Minimal Version (`minimal_predict.py`)

### Use When:
✅ Cell types DON'T match (e.g., human model → chimp data)  
✅ Want to find which cell types are functionally similar  
✅ Need to compare ALL predicted vs ALL true cell types  
✅ Don't care about matching cell type names

### What It Does:
```
Model (human)           Data (chimp)           Result
T_cell, B_cell    →     CD4_T, CD8_T, Mono →   Correlation matrix
   (3 types)              (3 types)            showing matches
                                               (3 × 3 matrix)
```

### Code:
```python
from minimal_predict import predict_accessibility, compare_across_celltypes, find_best_matches

adata_pred = predict_accessibility(model, adata)
corr = compare_across_celltypes(adata_pred, adata)
matches = find_best_matches(corr)
```

### Output:
- Predictions for MODEL's cell types (not input cell types)
- Correlation matrix: predicted × true cell types
- Best matches: which cell types correspond

---

## Full Version (`crested_utils.py`)

### Use When:
✅ Cell types DO match (or can be aligned)  
✅ Want predictions for SAME cell types as input  
✅ Need consistent dimensions for comparison  
✅ Have large datasets (need chunking)

### What It Does:
```
Model (human)           Data (chimp)           Result
T_cell, B_cell    →     T_cell, B_cell  →      Same cell types
   (2 types)              (2 types)              (2 types)
                          + align missing         + predictions
```

### Code:
```python
from crested_utils import predict_regions, compare_predictions

adata_pred = predict_regions(
    model, adata,
    target_cell_types=reference_celltypes,
    align_cell_types=True
)
corr_matrix, self_corr = compare_predictions(adata_pred)
```

### Output:
- Predictions for INPUT's cell types (aligned to reference)
- Correlation: predicted[i] vs true[i] (diagonal)
- Same dimensions as input (or target)

---

## Side-by-Side Comparison

| Feature | Minimal | Full |
|---------|---------|------|
| **Cell type matching** | Not required | Required (or aligned) |
| **Output cell types** | Model's types | Input's types |
| **Comparison** | All pairs | Same cell type |
| **Use case** | Cross-species | Same species |
| **Alignment** | Not needed | Automatic |
| **Correlation** | n_model × n_input | n_input × n_input |
| **Memory** | Basic | Optimized (chunking) |
| **Code lines** | ~300 | ~700 |
| **Functions** | 3 | 6 |

---

## Visual Examples

### Minimal Version Workflow

```
Step 1: Predict
┌──────────────────────┐         ┌──────────────────────┐
│ Input: Chimp Data    │         │ Output: Predictions  │
│                      │         │                      │
│ Cell types:          │   ───►  │ Cell types:          │
│  - CD4_T             │         │  - T_cell (model)    │
│  - CD8_T             │         │  - B_cell (model)    │
│  - Monocyte          │         │  - NK_cell (model)   │
│                      │         │                      │
│ Regions: 1000        │         │ Regions: 1000        │
└──────────────────────┘         └──────────────────────┘
   Different types!                Model's types!

Step 2: Compare All Pairs
                 CD4_T   CD8_T   Monocyte
T_cell (model)   0.85    0.82    0.45
B_cell (model)   0.23    0.21    0.78
NK_cell (model)  0.67    0.71    0.38

Step 3: Find Matches
T_cell matches CD4_T (r=0.85)
B_cell matches Monocyte (r=0.78)
NK_cell matches CD8_T (r=0.71)
```

### Full Version Workflow

```
Step 1: Align
┌──────────────────────┐         ┌──────────────────────┐
│ Input: Chimp Data    │         │ Target: Human Types  │
│                      │         │                      │
│ Cell types:          │         │ Cell types:          │
│  - T_cell            │   ───►  │  - T_cell            │
│  - B_cell            │         │  - B_cell            │
│  (missing NK_cell)   │         │  - NK_cell (→ zeros) │
│                      │         │                      │
│ Regions: 1000        │         │ Regions: 1000        │
└──────────────────────┘         └──────────────────────┘
   Same names!                     Aligned!

Step 2: Predict
Input cell types = Output cell types
(aligned to target)

Step 3: Compare Same Cell Types
                 Predicted
True    T_cell   B_cell   NK_cell
T_cell  0.92     0.15     0.23
B_cell  0.18     0.88     0.12
NK_cell 0.25     0.11     0.00 (missing)

Step 4: Self-Correlation
T_cell:  0.92 ✓ Good!
B_cell:  0.88 ✓ Good!
NK_cell: 0.00 ✗ Missing (filled with zeros)
```

---

## Real-World Scenarios

### Scenario 1: Cross-Species (Zero Overlap)
**Use: MINIMAL**

```python
# Human model, chimp data - cell types don't match
from minimal_predict import *

adata_pred = predict_accessibility(human_model, chimp_adata)
corr = compare_across_celltypes(adata_pred, chimp_adata)
matches = find_best_matches(corr)

# Result: Human T_cell matches Chimp CD4_T (r=0.85)
```

### Scenario 2: Same Species, Different Batch
**Use: FULL**

```python
# Both human, same cell type names
from crested_utils import predict_regions, compare_predictions

adata_pred = predict_regions(model, adata_batch2)
corr_matrix, self_corr = compare_predictions(adata_pred)

# Result: T_cell predicted with r=0.92 correlation
```

### Scenario 3: Cross-Species but Mapped Cell Types
**Use: FULL** (with alignment)

```python
# Human model, chimp data - but you know the mapping
from crested_utils import predict_regions

# Define mapping
chimp_to_human = {
    'CD4_T': 'T_cell',
    'CD8_T': 'T_cell',
    'Mono': 'Macrophage',
    ...
}

# Rename chimp cell types to human
adata_chimp.obs_names = [chimp_to_human[ct] for ct in adata_chimp.obs_names]

# Now use full version with alignment
adata_pred = predict_regions(
    model, adata_chimp,
    target_cell_types=human_celltypes,
    align_cell_types=True
)
```

### Scenario 4: Discover Cell Type Relationships
**Use: MINIMAL**

```python
# You DON'T know which cell types correspond
# Let the correlation matrix tell you!
from minimal_predict import *

adata_pred = predict_accessibility(model_speciesA, adata_speciesB)
corr = compare_across_celltypes(adata_pred, adata_speciesB)

# Heatmap shows which cell types are similar
import seaborn as sns
sns.heatmap(corr, annot=True)
```

---

## Quick Decision Guide

**Start here:**

```
Q: Do you know which cell types correspond?
│
├─ NO → Use MINIMAL
│   └─ Let correlation find matches
│
└─ YES → Can you rename/map them?
    │
    ├─ NO → Use MINIMAL
    │   └─ Names too different to map
    │
    └─ YES → Use FULL with alignment
        └─ Rename or use target_cell_types
```

---

## Code Comparison

### Your Goal: Predict and Compare

**Minimal (3 lines):**
```python
adata_pred = predict_accessibility(model, adata)
corr = compare_across_celltypes(adata_pred, adata)
matches = find_best_matches(corr)
```

**Full (2 lines):**
```python
adata_pred = predict_regions(model, regions=adata, target_cell_types=ref_types, align_cell_types=True)
corr_matrix, self_corr = compare_predictions(adata_pred)
```

---

## Which One for You?

Based on your original question: **"zero overlap in cell types"**

→ **Use MINIMAL version** (`minimal_predict.py`)

```python
from minimal_predict import predict_accessibility, compare_across_celltypes, find_best_matches

# Model trained on species A
# Data from species B with DIFFERENT cell types
adata_pred = predict_accessibility(model_A, adata_B)

# Compare ALL predicted vs ALL true
corr = compare_across_celltypes(adata_pred, adata_B)

# Find which cell types match
matches = find_best_matches(corr)
print(matches)
```

This will:
1. ✅ Work with zero cell type overlap
2. ✅ Show which cell types are functionally similar
3. ✅ Give correlation for all pairs
4. ✅ Help identify cell type correspondences

---

## Summary

| Your Situation | Use This |
|----------------|----------|
| Cell types don't match at all | **MINIMAL** |
| Cross-species, different annotations | **MINIMAL** |
| Want to discover cell type matches | **MINIMAL** |
| Cell types match by name | **FULL** |
| Same species/batch | **FULL** |
| Can map cell types manually | **FULL** (with rename) |
| Very large dataset (>10k regions) | **FULL** (has chunking) |

**For your specific case (zero overlap):** Use **`minimal_predict.py`** ✨
