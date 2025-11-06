"""
CREsted utility functions for cross-species enhancer predictions.

This package provides tools for predicting genomic accessibility and comparing
predictions with real scATAC-seq data using the CREsted framework.
"""

from .crested_utils import (
    predict_regions,
    align_adata_cell_types,
    compare_predictions,
    rowwise_correlation,
    resize_region,
)

from .minimal_predict import (
    predict_accessibility,
    compare_across_celltypes,
    find_best_matches,
)

from .insilico_mutagenesis_vect import (
    insilico_mutagenesis_vect,
    snp_mutagenesis_from_bed,
)

__all__ = [
    # Full version (crested_utils.py)
    "predict_regions",
    "align_adata_cell_types",
    "compare_predictions",
    "rowwise_correlation",
    "resize_region",
    # Minimal version (minimal_predict.py)
    "predict_accessibility",
    "compare_across_celltypes",
    "find_best_matches",
    # In-silico mutagenesis
    "insilico_mutagenesis_vect",
    "snp_mutagenesis_from_bed",
]

__version__ = "1.0.0"
