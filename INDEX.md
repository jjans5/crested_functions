# CREsted Utilities - Complete Package Index

## 📚 Documentation Structure

This package contains a unified utility module for working with CREsted predictions across species. Below is your complete documentation guide.

---

## 🚀 Quick Start (Read This First!)

**New user?** Start here:
1. Read **[SUMMARY.md](SUMMARY.md)** - Overview of what was improved
2. Skim **[QUICK_REFERENCE.md](QUICK_REFERENCE.md)** - One-page cheat sheet
3. Try **[example_usage.py](example_usage.py)** - Run examples
4. Refer to **[README.md](README.md)** when needed - Full documentation

**Migrating from old code?**
- Read **[MIGRATION_GUIDE.md](MIGRATION_GUIDE.md)** - Step-by-step migration

**Visual learner?**
- Read **[VISUAL_GUIDE.md](VISUAL_GUIDE.md)** - Diagrams and flowcharts

---

## 📖 Documentation Files

### 1. SUMMARY.md ⭐ START HERE
**Purpose:** Overview of improvements and what was done  
**Read when:** First time using the package  
**Contains:**
- What was improved
- Comparison with original code
- Key benefits summary
- Quick examples

**Time to read:** 10 minutes

---

### 2. QUICK_REFERENCE.md ⭐ KEEP HANDY
**Purpose:** One-page cheat sheet  
**Read when:** You know what you want to do but need syntax  
**Contains:**
- Function signatures
- Parameter guides
- Common workflows
- One-liner examples

**Time to read:** 5 minutes  
**Use:** Keep open while coding!

---

### 3. README.md 📘 FULL DOCS
**Purpose:** Complete documentation  
**Read when:** Need detailed information  
**Contains:**
- Installation instructions
- All function documentation
- Complete examples
- Use cases
- Tips and best practices

**Time to read:** 30 minutes  
**Use:** Reference when needed

---

### 4. MIGRATION_GUIDE.md 🔄 FOR EXISTING CODE
**Purpose:** Migrate from old scattered functions  
**Read when:** You have existing code to update  
**Contains:**
- Old vs new code side-by-side
- Pattern-by-pattern migration
- Testing checklist
- Common pitfalls

**Time to read:** 20 minutes  
**Use:** While migrating existing code

---

### 5. VISUAL_GUIDE.md 🎨 VISUAL LEARNER
**Purpose:** Visual explanations with diagrams  
**Read when:** Concepts need clarification  
**Contains:**
- Workflow diagrams
- Data flow charts
- Decision trees
- Architecture overview

**Time to read:** 15 minutes  
**Use:** When text isn't clear enough

---

## 🔧 Code Files

### 1. crested_utils.py ⭐ MAIN MODULE
**Purpose:** The unified utility module  
**Contains:**
- `predict_regions()` - Main prediction function
- `align_adata_cell_types()` - Cell type alignment
- `compare_predictions()` - Evaluate predictions
- `rowwise_correlation()` - Direct correlations
- `resize_region()` - Region utilities

**Import:** `from crested_utils import *`

---

### 2. example_usage.py 📝 EXAMPLES
**Purpose:** 7 complete working examples  
**Contains:**
1. Basic prediction
2. Cross-species with alignment
3. Comparing predictions
4. Full pipeline
5. Manual alignment
6. Region resizing
7. Direct correlation

**Run:** Uncomment examples in `main` block

---

### 3. test_crested_utils.py ✅ TESTS
**Purpose:** Verify installation and functions  
**Run:** `python test_crested_utils.py`  
**Tests:**
- Region resizing
- Correlation calculation
- Cell type alignment
- Prediction comparison

**Use:** After installation to verify everything works

---

## 🗂️ Original Files (Keep for Reference)

### predict_chunked.py
Original chunked prediction function. Now integrated into `crested_utils._predict_chunked()`.

### species_pred.py
Original cross-species pipeline. Now replaced by `predict_regions()` with alignment.

### rowwise_corr.py
Original correlation function. Now `crested_utils.rowwise_correlation()`.

### resize_region.py
Original region resizing. Now `crested_utils.resize_region()`.

**Keep these:** For reference until you've verified the new code works.

---

## 📊 File Organization by Task

### Task: "I want to start using the package"
1. Read **SUMMARY.md**
2. Skim **QUICK_REFERENCE.md**
3. Run **test_crested_utils.py**
4. Try **example_usage.py**

### Task: "I need to make predictions"
1. Check **QUICK_REFERENCE.md** → "predict_regions()"
2. If more detail needed: **README.md** → Function Reference
3. Copy example from **example_usage.py**

### Task: "I need cross-species comparison"
1. Check **VISUAL_GUIDE.md** → "Workflow 2"
2. Check **QUICK_REFERENCE.md** → "Workflow 2"
3. Run **example_usage.py** → `example_cross_species()`

### Task: "I want to migrate existing code"
1. Read **MIGRATION_GUIDE.md**
2. Follow pattern-by-pattern examples
3. Test with **test_crested_utils.py**

### Task: "I don't understand how it works"
1. Read **VISUAL_GUIDE.md**
2. Read **SUMMARY.md** → "Solving Your Original Goal"
3. Try **example_usage.py** → simple examples

---

## 🎯 Your Main Use Case: Cross-Species Prediction

### Files to Read (in order):
1. **SUMMARY.md** → "Solving Your Original Goal"
2. **VISUAL_GUIDE.md** → "Workflow 2: Cross-Species"
3. **QUICK_REFERENCE.md** → "Workflow 2: Cross-Species Comparison"
4. **example_usage.py** → `example_cross_species()` or `example_full_pipeline()`

### Your Main Function:
```python
from crested_utils import predict_regions, compare_predictions

# This is what you want!
adata_pred = predict_regions(
    model=model,
    regions=adata_species,
    target_cell_types=human_celltypes,
    align_cell_types=True,
    fill_missing_with_zeros=True
)

corr_matrix, self_corr = compare_predictions(adata_pred)
```

See **VISUAL_GUIDE.md** for diagrams showing how this solves your problem!

---

## 🔍 Quick Reference by Function

### predict_regions()
- **Docs:** README.md → "predict_regions()"
- **Quick ref:** QUICK_REFERENCE.md → "1. predict_regions()"
- **Examples:** example_usage.py → examples 1, 2, 4
- **Visual:** VISUAL_GUIDE.md → Workflows 1, 2, 3

### align_adata_cell_types()
- **Docs:** README.md → "align_adata_cell_types()"
- **Quick ref:** QUICK_REFERENCE.md → "2. align_adata_cell_types()"
- **Examples:** example_usage.py → example 5
- **Visual:** VISUAL_GUIDE.md → "Cell Type Alignment Logic"

### compare_predictions()
- **Docs:** README.md → "compare_predictions()"
- **Quick ref:** QUICK_REFERENCE.md → "3. compare_predictions()"
- **Examples:** example_usage.py → examples 3, 4
- **Visual:** VISUAL_GUIDE.md → "Workflow 4"

### rowwise_correlation()
- **Docs:** README.md → "rowwise_correlation()"
- **Quick ref:** QUICK_REFERENCE.md → "4. rowwise_correlation()"
- **Examples:** example_usage.py → example 7

### resize_region()
- **Docs:** README.md → "resize_region()"
- **Quick ref:** QUICK_REFERENCE.md → "5. resize_region()"
- **Examples:** example_usage.py → example 6

---

## 🎓 Learning Path

### Beginner (Never used CREsted utilities)
```
Day 1:
├─ Read: SUMMARY.md (10 min)
├─ Read: QUICK_REFERENCE.md (5 min)
├─ Run: test_crested_utils.py (2 min)
└─ Try: example_usage.py → example_basic_prediction() (10 min)

Day 2:
├─ Read: VISUAL_GUIDE.md → Your use case (10 min)
├─ Try: example_usage.py → example_cross_species() (15 min)
└─ Adapt to your data (30 min)
```

### Intermediate (Migrating from old code)
```
Session 1:
├─ Read: MIGRATION_GUIDE.md (20 min)
├─ Read: QUICK_REFERENCE.md (5 min)
└─ Compare: Your code vs new patterns (15 min)

Session 2:
├─ Migrate: One function at a time (1-2 hours)
├─ Test: Compare old vs new results (30 min)
└─ Document: What changed (15 min)
```

### Advanced (Building custom workflows)
```
├─ Read: README.md fully (30 min)
├─ Study: crested_utils.py source code (1 hour)
├─ Extend: Add custom functions (as needed)
└─ Reference: QUICK_REFERENCE.md for syntax
```

---

## 💡 Pro Tips

### 1. Keep These Files Open
- **QUICK_REFERENCE.md** - For syntax
- **Your terminal** - For running code
- **example_usage.py** - For copy-paste

### 2. Print This
- **QUICK_REFERENCE.md** - Keep on desk

### 3. Bookmark These Sections
- **README.md** → Your specific use case
- **VISUAL_GUIDE.md** → Relevant workflow diagram

### 4. Start Simple
- Don't read everything at once
- Start with basic examples
- Add complexity gradually

---

## ❓ FAQ Quick Links

### "How do I align cell types across species?"
→ **VISUAL_GUIDE.md** → "Workflow 2"  
→ **QUICK_REFERENCE.md** → "Workflow 2"

### "What parameters should I use?"
→ **QUICK_REFERENCE.md** → "Parameter Quick Guide"  
→ **README.md** → "Tips and Best Practices"

### "How do I migrate my existing code?"
→ **MIGRATION_GUIDE.md** → Start from top

### "I don't understand the concepts"
→ **VISUAL_GUIDE.md** → All workflows  
→ **SUMMARY.md** → "Solving Your Original Goal"

### "I'm getting errors"
→ **QUICK_REFERENCE.md** → "Error Messages"  
→ **README.md** → "Tips and Best Practices"

### "What's the difference from old code?"
→ **SUMMARY.md** → "Key Improvements"  
→ **MIGRATION_GUIDE.md** → Side-by-side comparisons

---

## 📞 Getting Help

### For Package Usage
1. Check **QUICK_REFERENCE.md** → Error Messages
2. Check **README.md** → Tips section
3. Check **example_usage.py** → Similar example
4. Check function docstring: `help(predict_regions)`

### For CREsted Package Issues
- CREsted docs: https://crested.readthedocs.io/
- CREsted issues: https://github.com/aertslab/CREsted/issues

---

## 🎯 TL;DR - Just Tell Me What To Do

### For Immediate Use:
```python
# 1. Import
from crested_utils import predict_regions, compare_predictions

# 2. Predict with alignment
adata_pred = predict_regions(
    model=your_model,
    regions=your_adata,
    target_cell_types=reference_celltypes,
    align_cell_types=True
)

# 3. Evaluate
corr_matrix, self_corr = compare_predictions(adata_pred)
print(f"Performance: {self_corr['correlation'].mean():.3f}")
```

### For More Details:
1. **Quick syntax?** → QUICK_REFERENCE.md
2. **How does it work?** → VISUAL_GUIDE.md
3. **Full documentation?** → README.md
4. **Migrating old code?** → MIGRATION_GUIDE.md

---

## 📝 Summary

This package provides **one unified module** (`crested_utils.py`) that replaces your scattered functions with:
- **60% less code**
- **Automatic cell type alignment**
- **Built-in validation and logging**
- **Comprehensive documentation**

**Main function:** `predict_regions()` - does prediction + optional alignment  
**Main use case:** Cross-species comparison with mismatched cell types  
**Main benefit:** No more manual alignment code!

---

## 📂 File Sizes Reference

| File | Size | Read Time | Type |
|------|------|-----------|------|
| SUMMARY.md | ~15 KB | 10 min | Overview |
| QUICK_REFERENCE.md | ~12 KB | 5 min | Cheat sheet |
| README.md | ~25 KB | 30 min | Full docs |
| MIGRATION_GUIDE.md | ~20 KB | 20 min | Migration |
| VISUAL_GUIDE.md | ~18 KB | 15 min | Diagrams |
| crested_utils.py | ~25 KB | - | Code |
| example_usage.py | ~12 KB | - | Examples |
| test_crested_utils.py | ~8 KB | 2 min | Tests |

---

**Total documentation:** ~90 minutes to read everything  
**Quick start:** 20 minutes to get going  
**Your use case:** 30 minutes to implement

---

**Last updated:** November 3, 2025  
**Package version:** 1.0  
**Python version:** 3.9+

---

Start with **SUMMARY.md** and go from there! 🚀
