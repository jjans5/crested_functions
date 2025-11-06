# Visual Workflow Guide

## Overview: What the Module Does

```
┌─────────────────────────────────────────────────────────────┐
│                     crested_utils.py                         │
│                                                              │
│  Input: AnnData or Regions → Predictions → Aligned Output   │
│                                                              │
│  Goal: Compare predictions across species with different    │
│        cell type compositions                                │
└─────────────────────────────────────────────────────────────┘
```

---

## Workflow 1: Basic Prediction

```
Input AnnData                predict_regions()              Output AnnData
┌──────────────┐                                          ┌──────────────┐
│ Cell Types:  │            ┌────────────┐               │ Cell Types:  │
│  - T_cell    │            │            │               │  - T_cell    │
│  - B_cell    │  ────────► │   Model    │  ───────────► │  - B_cell    │
│  - NK_cell   │            │            │               │  - NK_cell   │
│              │            └────────────┘               │              │
│ Regions:     │                                          │ Regions:     │
│  1000 peaks  │                                          │  1000 peaks  │
│              │                                          │              │
│ .X           │                                          │ .X           │
│ (true vals)  │                                          │ (true vals)  │
│              │                                          │ .layers['predicted']│
└──────────────┘                                          │ (predictions)│
                                                          └──────────────┘
```

**Code:**
```python
adata_pred = predict_regions(model, regions=adata)
```

---

## Workflow 2: Cross-Species with Alignment (YOUR MAIN USE CASE)

### Problem: Mismatched Cell Types

```
Human (Reference)              Chimpanzee (Species)
┌──────────────────┐          ┌──────────────────┐
│ Cell Types:      │          │ Cell Types:      │
│  1. T_cell       │          │  1. T_cell       │
│  2. B_cell       │          │  2. B_cell       │
│  3. NK_cell      │          │  3. NK_cell      │
│  4. Macrophage   │          │  4. Macrophage   │
│  5. Dendritic    │          │  5. Dendritic    │
│  6. Epithelial   │          │  (missing)       │  ← Not present!
│  7. Endothelial  │          │  (missing)       │  ← Not present!
└──────────────────┘          └──────────────────┘

Problem: Can't directly compare predictions!
         Dimensions don't match: 7 vs 5
```

### Solution: Automatic Alignment

```
Chimpanzee Data              predict_regions()             Aligned Output
┌──────────────────┐        with alignment              ┌──────────────────┐
│ Cell Types: 5    │                                    │ Cell Types: 7    │
│  - T_cell        │        ┌────────────┐             │  1. T_cell       │
│  - B_cell        │        │            │             │  2. B_cell       │
│  - NK_cell       │  ────► │   Model    │  ─────────► │  3. NK_cell      │
│  - Macrophage    │        │            │             │  4. Macrophage   │
│  - Dendritic     │        └────────────┘             │  5. Dendritic    │
│                  │              │                     │  6. Epithelial   │ ← Filled with 0
│ Regions: 1000    │              │                     │  7. Endothelial  │ ← Filled with 0
└──────────────────┘              │                     │                  │
                                  │                     │ Regions: 1000    │
                      target_cell_types                 └──────────────────┘
                    = [Human cell types]                
                                                        Now matches human!
                                                        Easy to compare!
```

**Code:**
```python
adata_pred = predict_regions(
    model=human_model,
    regions=adata_chimp,
    target_cell_types=human_celltypes,  # Reference
    align_cell_types=True,               # Do alignment
    fill_missing_with_zeros=True         # Fill missing
)
```

---

## Workflow 3: Multi-Species Comparison Pipeline

```
┌─────────────────────────────────────────────────────────────────────┐
│                     Multi-Species Analysis                           │
└─────────────────────────────────────────────────────────────────────┘

Step 1: Define Reference
┌──────────────────────────────────┐
│ Human Cell Types (Reference)     │
│  1. T_cell                       │
│  2. B_cell                       │
│  3. NK_cell                      │
│  4. Macrophage                   │
│  5. Dendritic                    │
│  6. Epithelial                   │
│  7. Endothelial                  │
└──────────────────────────────────┘
                │
                │ Apply to all species
                ▼
      ┌─────────┴─────────┐
      │                   │
┌─────▼──────┐    ┌──────▼──────┐    ┌──────▼──────┐
│ Chimpanzee │    │   Gorilla   │    │   Macaque   │
│            │    │             │    │             │
│ Has: 5/7   │    │ Has: 6/7    │    │ Has: 4/7    │
│ Missing: 2 │    │ Missing: 1  │    │ Missing: 3  │
└────────────┘    └─────────────┘    └─────────────┘
      │                   │                  │
      │ predict_regions() │                  │
      │ with alignment    │                  │
      ▼                   ▼                  ▼
┌────────────┐    ┌─────────────┐    ┌─────────────┐
│ Aligned    │    │ Aligned     │    │ Aligned     │
│ 7 cell     │    │ 7 cell      │    │ 7 cell      │
│ types      │    │ types       │    │ types       │
│ (missing   │    │ (missing    │    │ (missing    │
│  = zeros)  │    │  = zeros)   │    │  = zeros)   │
└────────────┘    └─────────────┘    └─────────────┘
      │                   │                  │
      │ compare_predictions()                │
      ▼                   ▼                  ▼
┌────────────┐    ┌─────────────┐    ┌─────────────┐
│ Corr: 0.85 │    │ Corr: 0.82  │    │ Corr: 0.78  │
└────────────┘    └─────────────┘    └─────────────┘

Now all species have same dimensions → Easy comparison!
```

**Code:**
```python
results = {}
for species in ['chimp', 'gorilla', 'macaque']:
    adata = load_species_data(species)
    
    # Predict with alignment
    adata_pred = predict_regions(
        model=human_model,
        regions=adata,
        target_cell_types=human_celltypes,
        align_cell_types=True
    )
    
    # Compare
    corr_matrix, self_corr = compare_predictions(adata_pred)
    results[species] = self_corr['correlation'].mean()

# All species now comparable!
for sp, corr in results.items():
    print(f"{sp}: {corr:.3f}")
```

---

## Workflow 4: Evaluation Pipeline

```
AnnData with Predictions         compare_predictions()         Results
┌─────────────────────┐                                    ┌──────────────────┐
│ .X (True values)    │                                    │ Correlation      │
│ ┌─────────────────┐ │         ┌──────────────┐         │ Matrix           │
│ │ T_cell    0.5...│ │         │              │         │ ┌──────────────┐ │
│ │ B_cell    0.3...│ │         │  Calculate   │         │ │ T→T  B→T ... │ │
│ │ NK_cell   0.2...│ │ ──────► │  Row-wise    │ ──────► │ │ T→B  B→B ... │ │
│ │ ...             │ │         │  Correlation │         │ │ ...          │ │
│ └─────────────────┘ │         │              │         │ └──────────────┘ │
│                     │         └──────────────┘         └──────────────────┘
│ .layers['predicted']│                 │                         │
│ ┌─────────────────┐ │                 │                         │
│ │ T_cell    0.4...│ │                 ▼                         ▼
│ │ B_cell    0.25..│ │         Self-Correlations        Statistics
│ │ NK_cell   0.18..│ │         ┌───────────────┐       ┌────────────────┐
│ │ ...             │ │         │ T_cell:  0.92 │       │ Mean:   0.85   │
│ └─────────────────┘ │         │ B_cell:  0.88 │       │ Median: 0.86   │
└─────────────────────┘         │ NK_cell: 0.81 │       │ Min:    0.75   │
                                 │ ...           │       │ Max:    0.95   │
                                 └───────────────┘       └────────────────┘
```

**Code:**
```python
corr_matrix, self_corr = compare_predictions(adata_pred)

print(f"Mean: {self_corr['correlation'].mean():.3f}")
print("\nPer cell type:")
print(self_corr.sort_values('correlation', ascending=False))

# Visualize
sns.heatmap(corr_matrix, cmap='RdBu_r', center=0)
```

---

## Data Flow Diagram

```
┌──────────────────────────────────────────────────────────────────┐
│                         YOUR WORKFLOW                             │
└──────────────────────────────────────────────────────────────────┘

Step 1: Load Data
┌───────────┐   ┌───────────┐   ┌───────────┐
│  Human    │   │ Chimp     │   │ Gorilla   │
│  AnnData  │   │ AnnData   │   │ AnnData   │
└─────┬─────┘   └─────┬─────┘   └─────┬─────┘
      │               │               │
      │ Train Model   │               │
      ▼               │               │
┌───────────┐         │               │
│  Trained  │         │               │
│   Model   │         │               │
└─────┬─────┘         │               │
      │               │               │
      │               │               │
      └───────┬───────┴───────┬───────┘
              │               │
              ▼               ▼

Step 2: Predict with Alignment
┌─────────────────────────────────────┐
│      predict_regions()               │
│  - Input: species AnnData            │
│  - target_cell_types: human types    │
│  - align_cell_types: True            │
│  - Output: aligned predictions       │
└─────────────────┬───────────────────┘
                  │
                  ▼

Step 3: Compare
┌─────────────────────────────────────┐
│      compare_predictions()           │
│  - Compare .X vs .layers['predicted']│
│  - Calculate correlations            │
│  - Output: matrix + statistics       │
└─────────────────┬───────────────────┘
                  │
                  ▼

Step 4: Visualize
┌─────────────────────────────────────┐
│         Heatmaps + Plots             │
│  - Correlation matrices              │
│  - Bar plots of self-correlations    │
│  - Cross-species comparisons         │
└─────────────────────────────────────┘
```

---

## Function Call Hierarchy

```
predict_regions()  ← YOUR MAIN ENTRY POINT
    │
    ├─► _predict_chunked()
    │       │
    │       ├─► crested.tl.predict()  [External]
    │       └─► ad.concat()            [External]
    │
    └─► align_adata_cell_types()  [If align=True]
            │
            └─► Creates aligned AnnData

compare_predictions()
    │
    └─► rowwise_correlation()
            │
            └─► Vectorized NumPy operations

resize_region()  [Standalone utility]
    │
    └─► String manipulation
```

---

## Memory Management Flow

```
Large Dataset (e.g., 50k regions)
         │
         ▼
┌─────────────────────────┐
│   Auto-detect memory    │
│   psutil.virtual_memory()│
└────────┬────────────────┘
         │
         ▼
   ┌────┴────┐
   │ Enough? │
   └────┬────┘
        │
   ┌────┴────────────┐
   │                 │
  Yes               No
   │                 │
   ▼                 ▼
Use large      Reduce chunk_size
chunks         (auto-adjust)
(fast)         (memory-safe)
   │                 │
   └────┬────────────┘
        │
        ▼
  Process in chunks
        │
        ├─► Chunk 1 (1000 regions)
        ├─► Chunk 2 (1000 regions)
        ├─► Chunk 3 (1000 regions)
        └─► ...
        │
        ▼
  Concatenate all chunks
        │
        ▼
  Return complete AnnData
```

---

## Cell Type Alignment Logic

```
Input Cell Types          Target Cell Types         Output
┌──────────────┐         ┌──────────────┐         ┌──────────────┐
│ T_cell       │ ─────┐  │ T_cell       │         │ T_cell       │ ✓ Present
│ B_cell       │ ──┐  │  │ B_cell       │         │ B_cell       │ ✓ Present
│ NK_cell      │ ┐ │  │  │ NK_cell      │         │ NK_cell      │ ✓ Present
│ Macrophage   │ │ │  │  │ Macrophage   │         │ Macrophage   │ ✓ Present
│              │ │ │  │  │ Dendritic    │         │ Dendritic    │ ✗ Missing → 0
│              │ │ │  │  │ Epithelial   │         │ Epithelial   │ ✗ Missing → 0
└──────────────┘ │ │  │  └──────────────┘         └──────────────┘
                 │ │  │         ▲
                 │ │  └─────────┤
                 │ └────────────┤
                 └──────────────┘
                        │
              align_adata_cell_types()
                fill_missing=True
```

---

## Quick Decision Tree

```
                    Start
                      │
                      ▼
        ┌─────────────────────────┐
        │ Need to make predictions?│
        └───────┬─────────────────┘
                │
          ┌─────┴─────┐
         Yes          No → Use other functions
          │
          ▼
    ┌───────────────────────┐
    │ Cross-species analysis?│
    └────┬──────────────────┘
         │
    ┌────┴────┐
   Yes       No
    │         │
    │         └─► predict_regions(model, adata)
    │
    ▼
┌─────────────────────────────────┐
│ predict_regions(                 │
│   model, adata,                  │
│   target_cell_types=ref_types,   │
│   align_cell_types=True          │
│ )                                │
└───────────┬─────────────────────┘
            │
            ▼
┌─────────────────────────────┐
│ Need to evaluate?            │
└─────┬───────────────────────┘
      │
  ┌───┴───┐
 Yes     No → Done!
  │
  ▼
┌──────────────────────────┐
│ compare_predictions(adata)│
└─────────┬────────────────┘
          │
          ▼
      Visualize!
```

---

## Key Concepts Summary

### 1. Chunking
```
Large dataset → Split into manageable chunks → Process each → Combine
                Memory-efficient!
```

### 2. Alignment
```
Different cell types → Add missing as zeros → Same dimensions → Comparable!
```

### 3. Row-wise Correlation
```
Cell type A (across regions) ↔ Cell type B (across regions)
How similar are their accessibility patterns?
```

### 4. Self-Correlation
```
True cell type A ↔ Predicted cell type A
How well did we predict this specific cell type?
```

---

## Visual Summary of Your Goal

```
BEFORE (Problem):
┌────────────────────────────────────────┐
│ Human:  7 cell types × 5000 regions    │
│ Chimp:  5 cell types × 4500 regions    │
│ Gorilla: 6 cell types × 4800 regions   │
│                                         │
│ ❌ Can't compare! Different dimensions  │
└────────────────────────────────────────┘

AFTER (Solution with crested_utils):
┌────────────────────────────────────────┐
│ Human:  7 cell types × common regions  │
│ Chimp:  7 cell types × common regions  │
│ Gorilla: 7 cell types × common regions │
│                                         │
│ ✅ Can compare! Same dimensions         │
│ ✅ Missing cell types filled with 0     │
│ ✅ Easy correlation calculation         │
└────────────────────────────────────────┘

Result:
┌────────────────────────────────────────┐
│ Chimp:  Mean correlation = 0.85        │
│ Gorilla: Mean correlation = 0.82       │
│ Macaque: Mean correlation = 0.78       │
│                                         │
│ ✅ Clear comparison across species!     │
└────────────────────────────────────────┘
```

---

This visual guide should help you understand how all the pieces fit together!
