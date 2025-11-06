# Getting Started with crested_utils

A 5-minute guide to get you up and running.

## Step 1: Installation (2 minutes)

```bash
# Install CREsted and dependencies
pip install crested anndata pandas numpy

# Optional: for memory monitoring
pip install psutil
```

## Step 2: Verify Installation (1 minute)

```bash
cd /path/to/crested_mod
python test_crested_utils.py
```

You should see:
```
🎉 All tests passed!
Your crested_utils module is working correctly.
```

## Step 3: Your First Prediction (2 minutes)

Create a file `my_analysis.py`:

```python
from crested_utils import predict_regions, compare_predictions
import crested
import anndata as ad

# Load your data and model
adata = ad.read_h5ad("your_data.h5ad")
model = crested.load_model("your_model")

# Make predictions
adata_pred = predict_regions(
    model=model,
    regions=adata,
    chunk_size=1000,
    verbose=True
)

# Evaluate
corr_matrix, self_corr = compare_predictions(adata_pred)
print(f"\nPrediction quality: {self_corr['correlation'].mean():.3f}")
print("\nPer cell type:")
print(self_corr.sort_values('correlation', ascending=False))
```

Run it:
```bash
python my_analysis.py
```

## Step 4: Cross-Species Analysis (5 minutes)

For your main use case (comparing across species):

```python
from crested_utils import predict_regions, compare_predictions
import crested
import anndata as ad

# Define reference cell types (e.g., from human)
human_celltypes = [
    'T_cell', 'B_cell', 'NK_cell', 'Macrophage', 
    'Dendritic_cell', 'Epithelial', 'Endothelial'
]

# Load model trained on human
model = crested.load_model("human_model.h5")

# Load species data (e.g., chimpanzee)
adata_chimp = ad.read_h5ad("chimpanzee_data.h5ad")

# Register genome
genome = crested.Genome("path/to/panTro3/genome.fa")
crested.register_genome(genome)

# Predict WITH automatic alignment
adata_pred = predict_regions(
    model=model,
    regions=adata_chimp,
    target_cell_types=human_celltypes,  # Align to human
    align_cell_types=True,               # Enable alignment
    fill_missing_with_zeros=True,        # Fill missing cell types
    chunk_size=2000,
    verbose=True
)

# Compare
corr_matrix, self_corr = compare_predictions(adata_pred)

print(f"\nChimp prediction quality: {self_corr['correlation'].mean():.3f}")

# Visualize
import matplotlib.pyplot as plt
import seaborn as sns

plt.figure(figsize=(10, 8))
sns.heatmap(corr_matrix, cmap='RdBu_r', center=0, vmin=-1, vmax=1, square=True)
plt.title(f"Chimp: Mean correlation = {self_corr['correlation'].mean():.3f}")
plt.tight_layout()
plt.savefig('chimp_correlation.pdf')
plt.show()
```

## That's It!

You're now using the improved utilities. 

### Next Steps

**If you want to learn more:**
- Read **QUICK_REFERENCE.md** for all functions
- Read **README.md** for complete documentation
- Check **example_usage.py** for 7 complete examples

**If you have existing code:**
- Read **MIGRATION_GUIDE.md** for how to migrate

**If concepts are unclear:**
- Read **VISUAL_GUIDE.md** for diagrams

**If you just need a reminder:**
- Keep **QUICK_REFERENCE.md** handy

### Common Issues

**"Import error: crested not found"**
```bash
pip install crested
```

**"Import error: crested_utils not found"**
Make sure you're in the correct directory or add it to your Python path:
```python
import sys
sys.path.append('/path/to/crested_mod')
from crested_utils import predict_regions
```

**"Out of memory errors"**
Reduce chunk_size:
```python
adata_pred = predict_regions(model, adata, chunk_size=500)
```

**"Cell types don't match"**
Use alignment:
```python
adata_pred = predict_regions(
    model, adata,
    target_cell_types=reference_celltypes,
    align_cell_types=True
)
```

### Quick Help

```python
# Get help on any function
help(predict_regions)
help(align_adata_cell_types)
help(compare_predictions)

# Or use your IDE's autocomplete with type hints!
```

### File Navigation

```
crested_mod/
├── crested_utils.py          ← Import from this
├── example_usage.py          ← Copy examples from this
├── test_crested_utils.py     ← Run this to test
├── QUICK_REFERENCE.md        ← Keep this open while coding
├── README.md                 ← Read when you need details
└── INDEX.md                  ← Navigate all documentation
```

---

## Your 5-Minute Checklist

- [ ] Install dependencies (`pip install crested anndata pandas`)
- [ ] Run tests (`python test_crested_utils.py`)
- [ ] Copy basic example from above
- [ ] Replace with your data/model paths
- [ ] Run it!
- [ ] Check QUICK_REFERENCE.md for more

---

**Need more help?** Start with **INDEX.md** to navigate all documentation.

**Ready to go?** You've got this! 🚀
