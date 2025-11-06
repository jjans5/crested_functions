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

from .mutation_scoring import (
    saturation_mutagenesis,
    score_snps,
    insilico_mutagenesis_vect,  # backward compatibility alias
    snp_mutagenesis_from_bed,   # backward compatibility alias
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
    # Mutation scoring (mutation_scoring.py)
    "saturation_mutagenesis",
    "score_snps",
    "insilico_mutagenesis_vect",  # backward compatibility
    "snp_mutagenesis_from_bed",   # backward compatibility
]

__version__ = "1.0.0"
