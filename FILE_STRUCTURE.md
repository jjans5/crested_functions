# Repository Structure

This document explains the organization of files in this repository.

## 📁 Core Modules

### Main Functionality
- **`crested_utils.py`** - Full-featured utilities with alignment and chunking
- **`minimal_predict.py`** - Minimal version for zero cell type overlap
- **`insilico_mutagenesis_vect.py`** - Vectorized in-silico mutagenesis

### Examples and Tests
- **`example_usage.py`** - 7 comprehensive examples for full version
- **`demo_minimal.py`** - Demo for minimal version with synthetic data
- **`insilico_mutagenesis_vect_example.py`** - ISM workflow example
- **`test_crested_utils.py`** - Automated tests

## 📚 Documentation

### Getting Started
- **`GITHUB_README.md`** - Main README for GitHub (rename to README.md when pushing)
- **`GETTING_STARTED.md`** - 5-minute quick start guide
- **`INDEX.md`** - Navigation guide for all documentation

### Reference Documentation
- **`QUICK_REFERENCE.md`** - One-page cheat sheet (keep handy!)
- **`README.md`** - Complete documentation for full version
- **`MINIMAL_README.md`** - Documentation for minimal version

### Guides
- **`WHICH_VERSION.md`** - Decision guide: minimal vs full
- **`VISUAL_GUIDE.md`** - Visual workflow diagrams
- **`MIGRATION_GUIDE.md`** - Migrate from older code
- **`SUMMARY.md`** - Overview of improvements

## 🗂️ Legacy Files (For Reference)

These are the original scattered functions that were consolidated:
- **`predict_chunked.py`** - Original chunked prediction
- **`species_pred.py`** - Original cross-species code
- **`rowwise_corr.py`** - Original correlation function
- **`resize_region.py`** - Region utilities

**Note:** These are kept for reference but you should use the new modules (`crested_utils.py` or `minimal_predict.py`) instead.

## 🔧 Configuration Files
- **`.gitignore`** - Git ignore patterns
- **`LICENSE`** - MIT License

## 📊 Recommended Reading Order

### For New Users
1. `GITHUB_README.md` (or `README.md` after push) - Overview
2. `GETTING_STARTED.md` - Quick start
3. `QUICK_REFERENCE.md` - Syntax reference
4. `example_usage.py` or `demo_minimal.py` - Try examples

### For Existing Users (Migrating)
1. `SUMMARY.md` - What changed
2. `MIGRATION_GUIDE.md` - How to migrate
3. `WHICH_VERSION.md` - Which version to use

### For Understanding Workflows
1. `VISUAL_GUIDE.md` - See diagrams
2. `WHICH_VERSION.md` - Decision trees
3. Examples: Try `demo_minimal.py` or `example_usage.py`

## 🎯 File Categories

### Essential (Must Read)
- `GITHUB_README.md` / `README.md`
- `GETTING_STARTED.md`
- `crested_utils.py` or `minimal_predict.py`

### Reference (Keep Handy)
- `QUICK_REFERENCE.md`
- `WHICH_VERSION.md`

### Detailed (Read When Needed)
- `README.md` (full docs)
- `MINIMAL_README.md`
- `VISUAL_GUIDE.md`
- `MIGRATION_GUIDE.md`

### Examples (Try Them)
- `demo_minimal.py`
- `example_usage.py`
- `insilico_mutagenesis_vect_example.py`

### Legacy (Reference Only)
- `predict_chunked.py`
- `species_pred.py`
- `rowwise_corr.py`
- `resize_region.py`

## 📦 What to Include in Your Project

### Minimal Setup
Just need the essentials:
```
your_project/
├── minimal_predict.py          # or crested_utils.py
├── QUICK_REFERENCE.md
└── GETTING_STARTED.md
```

### Full Setup
Complete toolkit:
```
your_project/
├── crested_utils.py
├── minimal_predict.py
├── insilico_mutagenesis_vect.py
├── QUICK_REFERENCE.md
├── README.md
├── MINIMAL_README.md
└── example_usage.py
```

### Development Setup
Everything including tests:
```
your_project/
├── crested_utils.py
├── minimal_predict.py
├── insilico_mutagenesis_vect.py
├── test_crested_utils.py
├── demo_minimal.py
├── All documentation files
└── Examples
```

## 🚀 For GitHub

When pushing to GitHub, the repository structure is ready to go. Just:

1. Rename `GITHUB_README.md` to `README.md` (or keep both)
2. Update your GitHub username in the README
3. Initialize git and push
4. The .gitignore is already configured

Users will see:
- Clear README on the main page
- Documentation in separate files
- Examples they can run
- Tests they can verify
- Original files for reference

---

This structure balances completeness with usability. Users can start quickly with the essentials, while having access to comprehensive documentation when needed.
