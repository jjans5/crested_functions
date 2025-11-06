# API Reference

## crested_utils

### predict_regions()

Predict accessibility for genomic regions.

**Parameters:**
- `model` - CREsted model
- `regions` - List of region strings OR AnnData object
- `genome` - Genome object (required if regions is a list)
- `batch_size` (default: 16) - Batch size for prediction
- `chunk_size` (default: 1000) - Regions per chunk
- `target_cell_types` - Cell types for alignment
- `align_cell_types` (default: True) - Align to target cell types
- `verbose` (default: True) - Enable logging

**Returns:** AnnData with predictions in `.layers['predicted']`

### align_adata_cell_types()

Align cell types between datasets.

**Parameters:**
- `adata` - Input AnnData
- `target_cell_types` - Desired cell types in order
- `fill_missing` (default: True) - Fill missing with zeros
- `verbose` (default: True) - Enable logging

**Returns:** AnnData with aligned cell types

### compare_predictions()

Compare predicted vs real values.

**Parameters:**
- `adata` - AnnData with predictions
- `prediction_layer` (default: 'predicted') - Layer name
- `true_layer` (default: None) - True values layer (None = use .X)
- `method` (default: 'pearson') - Correlation method

**Returns:** Dictionary with 'correlations' and 'diagonal_correlations'

### rowwise_correlation()

Compute row-wise correlations.

**Parameters:**
- `A` - First DataFrame
- `B` - Second DataFrame
- `method` (default: 'pearson') - Correlation method

**Returns:** Correlation matrix

## minimal_predict

### predict_accessibility()

Predict using model's own cell types (for zero overlap case).

**Parameters:**
- `model` - CREsted model
- `regions` - AnnData with regions
- `batch_size` (default: 16) - Batch size

**Returns:** AnnData with predictions

### compare_across_celltypes()

Compare all cell type pairs.

**Parameters:**
- `predicted_adata` - Predicted values
- `real_adata` - Real values
- `method` (default: 'pearson') - Correlation method

**Returns:** DataFrame with all pairwise correlations

### find_best_matches()

Find best matching cell types.

**Parameters:**
- `comparison_df` - Output from compare_across_celltypes()
- `top_n` (default: 5) - Number of matches per cell type

**Returns:** DataFrame with top matches

## insilico_mutagenesis_vect

### insilico_mutagenesis_vect()

Vectorized in-silico mutagenesis.

**Parameters:**
- `seq` - DNA sequence string
- `model` - CREsted model
- `adata` - AnnData with cell type names
- `chrom` - Chromosome name (optional)
- `start` - Genomic start position (optional)
- `batch_size` (default: 2) - Batch size
- `return_long` (default: True) - Return long format DataFrame

**Returns:** Dictionary with:
- `wt_pred` - Wildtype predictions (Series)
- `mut_df_wide` - All mutations in wide format (DataFrame)
- `mut_df_long` - Tidy format with log2fc and delta (DataFrame)

### snp_mutagenesis_from_bed()

Analyze SNPs from BED file.

**Parameters:**
- `bed_file` - Path to BED file
- `model` - CREsted model
- `adata` - AnnData with cell type names
- `genome` - Genome object
- `seq_length` - Sequence length (inferred from model if None)
- `batch_size` (default: 2) - Batch size
- `return_long` (default: True) - Return long format

**Returns:** DataFrame with columns:
- `chrom`, `pos`, `ref`, `alt`, `mut_id`
- `cell_type`, `wt_pred`, `mut_pred`
- `log2fc`, `delta`
