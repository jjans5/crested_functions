# CREsted Utilities

A collection of utility functions for working with the [CREsted](https://github.com/aertslab/CREsted) package, focusing on cross-species predictions and comparisons.

## 🚀 Overview

This repository provides enhanced utilities for:
- **Cross-species predictions**: Compare regulatory predictions across species with different cell types
- **In-silico mutagenesis**: Vectorized single-nucleotide mutagenesis for enhancer analysis
- **Cell type alignment**: Automatically align cell types when comparing across datasets
- **Accessibility comparison**: Compare predicted vs observed chromatin accessibility

## 📦 What's Included

### Two Versions Available

#### 1. **Minimal Version** (`minimal_predict.py`) - For Zero Cell Type Overlap
Perfect when model and data have completely different cell types (e.g., cross-species):

```python
from minimal_predict import predict_accessibility, compare_across_celltypes, find_best_matches

# Predict using model's cell types
adata_pred = predict_accessibility(model, adata)

# Compare all predicted vs all true cell types
corr = compare_across_celltypes(adata_pred, adata)

# Find functional matches
matches = find_best_matches(corr)
```

**Use when**: Cell types don't match (e.g., human model → chimp data)

#### 2. **Full Version** (`crested_utils.py`) - For Aligned Cell Types
Comprehensive utilities with automatic alignment and memory optimization:

```python
from crested_utils import predict_regions, compare_predictions

# Predict with automatic alignment
adata_pred = predict_regions(
    model=model,
    regions=adata,
    target_cell_types=reference_celltypes,
    align_cell_types=True
)

# Compare predictions
corr_matrix, self_corr = compare_predictions(adata_pred)
```

**Use when**: Cell types match or can be aligned

### In-Silico Mutagenesis

Vectorized single-nucleotide mutagenesis for efficient enhancer analysis:

```python
from insilico_mutagenesis_vect import insilico_mutagenesis_vect

result = insilico_mutagenesis_vect(
    seq=sequence,
    model=model,
    adata=adata,
    chrom="chr1",
    start=12345,
    return_long=True
)

# Get wildtype predictions
wt_pred = result["wt_pred"]

# Get all mutant predictions
mut_df_wide = result["mut_df_wide"]
mut_df_long = result["mut_df_long"]
```

## 📚 Documentation

- **[GETTING_STARTED.md](GETTING_STARTED.md)** - Quick 5-minute start guide
- **[QUICK_REFERENCE.md](QUICK_REFERENCE.md)** - One-page cheat sheet
- **[README.md](README.md)** - Complete documentation (full version)
- **[MINIMAL_README.md](MINIMAL_README.md)** - Documentation for minimal version
- **[WHICH_VERSION.md](WHICH_VERSION.md)** - Decision guide: which version to use
- **[VISUAL_GUIDE.md](VISUAL_GUIDE.md)** - Visual workflow diagrams
- **[MIGRATION_GUIDE.md](MIGRATION_GUIDE.md)** - Migrate from older code
- **[INDEX.md](INDEX.md)** - Navigate all documentation

## 🎯 Quick Start

### Installation

```bash
# Install CREsted and dependencies
pip install crested anndata pandas numpy

# Optional: memory monitoring
pip install psutil

# Clone this repository
git clone https://github.com/jjans5/crested-utils.git
cd crested-utils
```

### Basic Usage

#### Cross-Species Comparison (Zero Cell Type Overlap)

```python
from minimal_predict import predict_accessibility, compare_across_celltypes, find_best_matches
import crested
import anndata as ad

# Load model trained on human
model = crested.load_model("human_model.h5")

# Load chimp data (different cell types)
adata_chimp = ad.read_h5ad("chimp_data.h5ad")

# Predict
adata_pred = predict_accessibility(model, adata_chimp)

# Compare
corr = compare_across_celltypes(adata_pred, adata_chimp)

# Find matches
matches = find_best_matches(corr, top_k=3)
print(matches)
```

#### Same Species with Alignment

```python
from crested_utils import predict_regions, compare_predictions

# Predict with alignment
adata_pred = predict_regions(
    model=model,
    regions=adata,
    target_cell_types=['T_cell', 'B_cell', 'NK_cell'],
    align_cell_types=True
)

# Evaluate
corr_matrix, self_corr = compare_predictions(adata_pred)
print(f"Mean correlation: {self_corr['correlation'].mean():.3f}")
```

#### In-Silico Mutagenesis

```python
from insilico_mutagenesis_vect import insilico_mutagenesis_vect
import crested

# Load model and genome
model = crested.load_model("model.h5")
genome = crested.Genome("genome.fa")

# Fetch sequence
seq = genome.fetch("chr1", 100000, 102114)

# Run ISM
result = insilico_mutagenesis_vect(
    seq=seq,
    model=model,
    adata=adata,
    chrom="chr1",
    start=100000,
    return_long=True
)

# Analyze results
wt_pred = result["wt_pred"]
mut_df = result["mut_df_long"]

# Find most impactful mutations
top_impact = mut_df.sort_values("delta").head(20)
```

## 🧪 Testing

```bash
# Test basic functions
python test_crested_utils.py

# Try demo with synthetic data
python demo_minimal.py
```

## 📖 Examples

See [`example_usage.py`](example_usage.py) for 7 comprehensive examples:
1. Basic prediction
2. Cross-species with alignment
3. Comparing predictions
4. Full multi-species pipeline
5. Manual cell type alignment
6. Region resizing
7. Direct correlation analysis

See [`insilico_mutagenesis_vect_example.py`](insilico_mutagenesis_vect_example.py) for ISM workflow.

## 🔑 Key Features

### Minimal Version
✅ Works with **zero cell type overlap**  
✅ Finds **functional cell type matches**  
✅ Compares **all pairs** of cell types  
✅ **Ultra-simple**: 3 functions, 3 lines of code

### Full Version
✅ **Automatic cell type alignment**  
✅ **Memory-optimized** chunking  
✅ **Progress tracking** with detailed logging  
✅ **Type hints** for IDE support  
✅ **Comprehensive validation**

### In-Silico Mutagenesis
✅ **Vectorized** for speed  
✅ Returns both **wide and long** formats  
✅ Automatic **delta calculation**  
✅ **Genomic coordinates** support

## 📊 Use Cases

### Cross-Species Regulatory Analysis
Compare chromatin accessibility predictions between species where cell type annotations differ:
```python
# Human model, chimp/gorilla/macaque data
for species in ['chimp', 'gorilla', 'macaque']:
    adata = load_species_data(species)
    adata_pred = predict_accessibility(human_model, adata)
    corr = compare_across_celltypes(adata_pred, adata)
    # Identify functionally similar cell types
```

### Enhancer Variant Analysis
Systematically test all single-nucleotide variants in regulatory regions:
```python
# Test all SNPs in an enhancer
result = insilico_mutagenesis_vect(enhancer_seq, model, adata)
# Find mutations with largest impact
top_variants = result["mut_df_long"].sort_values("delta")
```

### Model Evaluation
Validate model predictions across different batches or conditions:
```python
# Compare predictions vs reality
adata_pred = predict_regions(model, test_data, align_cell_types=True)
corr_matrix, self_corr = compare_predictions(adata_pred)
```

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- Built for the [CREsted](https://github.com/aertslab/CREsted) package
- Developed for cross-species regulatory genomics analysis

## 📞 Citation

If you use these utilities in your research, please cite:

```bibtex
@software{crested_utils,
  author = {Janssens, J.},
  title = {CREsted Utilities: Tools for Cross-Species Regulatory Predictions},
  year = {2025},
  url = {https://github.com/jjans5/crested-utils}
}
```

And please also cite the CREsted package:

> Kempynck, N., De Winter, S., et al. CREsted: modeling genomic and synthetic cell type-specific enhancers across tissues and species. bioRxiv (2025).

## 🔗 Links

- [CREsted Documentation](https://crested.readthedocs.io/)
- [CREsted GitHub](https://github.com/aertslab/CREsted)
- [Issue Tracker](https://github.com/jjans5/crested-utils/issues)

## 📝 Quick Reference

| Task | Use This | Code |
|------|----------|------|
| Cross-species (no overlap) | Minimal | `predict_accessibility(model, adata)` |
| Same species | Full | `predict_regions(model, adata, align=True)` |
| In-silico mutagenesis | ISM | `insilico_mutagenesis_vect(seq, model, adata)` |
| Compare predictions | Full/Minimal | `compare_predictions()` / `compare_across_celltypes()` |

---

**Need help?** Start with [GETTING_STARTED.md](GETTING_STARTED.md) or check [INDEX.md](INDEX.md) to navigate all documentation.
