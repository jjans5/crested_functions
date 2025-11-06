# Minimal Predict - Zero Cell Type Overlap Solution

**Ultra-minimal code for predicting and comparing accessibility when model and data have ZERO overlapping cell types.**

## The Problem

You have:
- Model trained on **Species A** (e.g., human: T_cell, B_cell, NK_cell)
- Data from **Species B** (e.g., chimp: CD4_T, CD8_T, Monocyte)
- **Zero overlap** in cell type names!

How do you compare predictions vs reality?

## The Solution

**3 functions, 3 lines of code:**

```python
from minimal_predict import predict_accessibility, compare_across_celltypes, find_best_matches

# 1. Predict (model predicts for ITS cell types)
adata_pred = predict_accessibility(model, adata)

# 2. Compare (all predicted vs all true cell types)
corr_matrix = compare_across_celltypes(adata_pred, adata)

# 3. Analyze (find which cell types match)
matches = find_best_matches(corr_matrix)
```

**That's it!**

---

## How It Works

### Step 1: Predict

```python
adata_pred = predict_accessibility(model, adata_chimp)
```

**Input:** Chimp data with chimp cell types  
**Output:** Predictions for MODEL's cell types (e.g., human)

```
Chimp Data                Model (trained on human)        Predictions
┌──────────────┐         ┌──────────────┐               ┌──────────────┐
│ CD4_T        │         │              │               │ T_cell       │
│ CD8_T        │  ────►  │  Human Model │  ──────────►  │ B_cell       │
│ Monocyte     │         │              │               │ NK_cell      │
└──────────────┘         └──────────────┘               └──────────────┘
Chimp celltypes          Trained on human               Human celltypes
× chimp regions                                         × chimp regions
```

### Step 2: Compare

```python
corr_matrix = compare_across_celltypes(adata_pred, adata_chimp)
```

Computes correlation between **every** predicted cell type and **every** true cell type.

```
Correlation Matrix:
                 CD4_T   CD8_T   Monocyte   ...
T_cell           0.85    0.82    0.45       ...
B_cell           0.23    0.21    0.15       ...
NK_cell          0.67    0.71    0.38       ...
```

High correlation = similar accessibility patterns!

### Step 3: Analyze

```python
matches = find_best_matches(corr_matrix, top_k=3)
```

Finds best matching true cell type for each predicted cell type.

```
predicted_celltype   true_celltype   correlation   rank
T_cell               CD4_T           0.85          1
T_cell               CD8_T           0.82          2
T_cell               NK_cell_chimp   0.67          3
B_cell               B_cell_chimp    0.78          1
...
```

---

## Complete Example

```python
import crested
import anndata as ad
from minimal_predict import (
    predict_accessibility,
    compare_across_celltypes,
    find_best_matches
)

# Load model (trained on human)
model = crested.load_model("human_model.h5")

# Load data (chimp with different cell types)
adata_chimp = ad.read_h5ad("chimp_data.h5ad")

# Register genome
genome = crested.Genome("chimp_genome.fa")
crested.register_genome(genome)

# 1. Predict
adata_pred = predict_accessibility(model, adata_chimp)

# 2. Compare
corr_matrix = compare_across_celltypes(adata_pred, adata_chimp)

# 3. Analyze
matches = find_best_matches(corr_matrix, top_k=3)
print(matches)

# 4. Visualize
import matplotlib.pyplot as plt
import seaborn as sns

plt.figure(figsize=(10, 8))
sns.heatmap(corr_matrix, cmap='RdBu_r', center=0, annot=True, fmt='.2f')
plt.title('Model Cell Types vs True Cell Types')
plt.savefig('celltype_matching.pdf')
```

---

## Why This Works

**Key insight:** You don't need cell types to match by name!

- Model predicts accessibility **patterns** for its cell types
- You compare these patterns to your actual data
- High correlation = functionally similar cell types
- Even if names are completely different!

**Example:**
- Model's "T_cell" has high correlation with your "CD4_T" and "CD8_T"
- → They're functionally similar even though names differ!
- → Model's predictions are valid for your T cells

---

## Demo

Try the demo with synthetic data:

```bash
python demo_minimal.py
```

This shows the complete workflow without needing real data files.

---

## Functions

### `predict_accessibility(model, adata, batch_size=16)`

Predict accessibility for regions in adata.

**Returns:** AnnData with predictions for MODEL's cell types (not input cell types!)

### `compare_across_celltypes(adata_pred, adata_true, method='pearson')`

Compare predicted vs true across all cell type pairs.

**Returns:** Correlation matrix (predicted × true cell types)

### `find_best_matches(corr_matrix, top_k=3)`

Find best matching true cell types for each predicted cell type.

**Returns:** DataFrame with matches sorted by correlation

---

## Key Points

✅ **Works with zero cell type overlap**  
✅ **Model predicts for its own cell types**  
✅ **Compares patterns, not names**  
✅ **Finds functional matches**  
✅ **3 lines of code**

❌ Does NOT require cell types to match  
❌ Does NOT need manual alignment  
❌ Does NOT need same number of cell types

---

## Comparison with Full Package

**Minimal version** (`minimal_predict.py`):
- 3 functions, ~300 lines
- For when cell types DON'T match
- Compares ALL pairs
- Finds functional matches

**Full version** (`crested_utils.py`):
- 6 functions, ~700 lines
- For when cell types DO match (or can be aligned)
- Aligns cell types automatically
- More features (chunking, memory management, etc.)

**Use minimal version when:** Cell types have zero overlap (cross-species, different annotations, etc.)

**Use full version when:** Cell types match or can be aligned (same species, consistent annotations)

---

## Your Specific Use Case

**Scenario:** Human model → Chimp data

```python
# Load human model
human_model = crested.load_model("human_model.h5")

# Load chimp data (different cell types!)
adata_chimp = ad.read_h5ad("chimp_data.h5ad")

# Predict (model uses human cell types)
adata_pred = predict_accessibility(human_model, adata_chimp)
# Result: human cell types × chimp regions

# Compare (all human vs all chimp cell types)
corr = compare_across_celltypes(adata_pred, adata_chimp)
# Result: matrix showing which cell types match

# Find matches
matches = find_best_matches(corr)
# Result: "human T_cell matches chimp CD4_T with r=0.85"
```

**Interpretation:**
- High correlation = model's predictions are good for that cell type
- You can identify which human cell type corresponds to which chimp cell type
- Even though names are completely different!

---

## Files

- **`minimal_predict.py`** - Main module (3 functions)
- **`demo_minimal.py`** - Demo with synthetic data
- **`MINIMAL_README.md`** - This file

---

## Quick Start

```python
from minimal_predict import predict_accessibility, compare_across_celltypes, find_best_matches

adata_pred = predict_accessibility(model, adata)
corr = compare_across_celltypes(adata_pred, adata)
matches = find_best_matches(corr)
```

**Done!** 🎉
