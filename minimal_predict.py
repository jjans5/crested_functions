"""
Minimal code for predicting accessibility and comparing across different cell types.

Use case: Model trained on species A, predict on species B with different cell types.
Goal: Compare predicted vs true accessibility even with zero cell type overlap.
"""

import numpy as np
import pandas as pd
import anndata as ad
import crested
from typing import Optional, List
import warnings


def predict_accessibility(
    model,
    adata: ad.AnnData,
    batch_size: int = 16,
    layer_name: str = "predicted"
) -> ad.AnnData:
    """
    Predict accessibility for regions in adata using model.
    
    Works even if model was trained on completely different cell types.
    Model will predict for ITS OWN cell types, not the ones in adata.
    
    Parameters
    ----------
    model : CREsted model
        Trained model (e.g., on human cell types)
    adata : AnnData
        Data with regions to predict (e.g., chimp cell types)
    batch_size : int
        Batch size for prediction
    layer_name : str
        Name for output layer
        
    Returns
    -------
    AnnData
        New AnnData with:
        - .obs: MODEL's cell types (what model was trained on)
        - .var: Input regions (from adata)
        - .X: Predictions for model's cell types
        - .layers[layer_name]: Same as .X
        
    Example
    -------
    >>> # Model trained on human: T_cell, B_cell, NK_cell
    >>> # Data has chimp: CD4_T, CD8_T, B_cell, Monocyte
    >>> adata_pred = predict_accessibility(human_model, chimp_adata)
    >>> # Result has human cell types (T_cell, B_cell, NK_cell) 
    >>> # × chimp regions
    """
    
    # Get regions from input
    regions = adata.var_names.tolist()
    n_regions = len(regions)
    
    print(f"Input: {adata.n_obs} cell types × {n_regions} regions")
    print(f"Predicting accessibility using model...")
    
    # Create minimal adata for prediction
    # We only need the regions - model will predict for its own cell types
    minimal_adata = ad.AnnData(
        X=np.zeros((1, n_regions), dtype=np.float32),  # Dummy data
        var=adata.var.copy()
    )
    minimal_adata.obs['dummy'] = ['placeholder']
    
    # Predict - this returns (n_regions, n_model_celltypes)
    predictions = crested.tl.predict(
        minimal_adata,
        model,
        batch_size=batch_size
    )
    
    # predictions shape: (n_regions, n_model_celltypes)
    # We want: (n_model_celltypes, n_regions)
    predictions = predictions.T
    
    # Get model's cell types from prediction shape
    n_model_celltypes = predictions.shape[0]
    model_celltype_names = [f"model_celltype_{i}" for i in range(n_model_celltypes)]
    
    print(f"Output: {n_model_celltypes} model cell types × {n_regions} regions")
    print(f"Note: Model predicts for its own cell types, not input cell types")
    
    # Create output AnnData
    obs = pd.DataFrame(index=model_celltype_names)
    var = adata.var.copy()
    
    adata_pred = ad.AnnData(X=predictions, obs=obs, var=var)
    adata_pred.layers[layer_name] = predictions
    
    # Store metadata about input
    adata_pred.uns['input_cell_types'] = adata.obs_names.tolist()
    adata_pred.uns['input_shape'] = adata.shape
    
    return adata_pred


def compare_across_celltypes(
    adata_pred: ad.AnnData,
    adata_true: ad.AnnData,
    method: str = "pearson"
) -> pd.DataFrame:
    """
    Compare predicted vs true accessibility across ALL cell type pairs.
    
    This computes correlation between every predicted cell type and every
    true cell type. Useful when cell types don't match (e.g., cross-species).
    
    Parameters
    ----------
    adata_pred : AnnData
        Predictions (model's cell types × regions)
    adata_true : AnnData
        True data (actual cell types × regions)
    method : str
        'pearson' or 'spearman'
        
    Returns
    -------
    DataFrame
        Correlation matrix (pred_celltypes × true_celltypes)
        Each value is correlation across common regions
        
    Example
    -------
    >>> # Model: human (T, B, NK) × chimp regions
    >>> # True: chimp (CD4_T, CD8_T, B, Mono) × chimp regions
    >>> corr = compare_across_celltypes(pred, true)
    >>> # Result: 3 × 4 matrix showing which cell types match
    """
    
    # Get common regions
    common_regions = adata_pred.var_names.intersection(adata_true.var_names)
    if len(common_regions) == 0:
        raise ValueError("No common regions between predicted and true data!")
    
    print(f"Comparing across {len(common_regions)} common regions")
    print(f"Predicted: {adata_pred.n_obs} cell types")
    print(f"True: {adata_true.n_obs} cell types")
    
    # Subset to common regions
    pred_subset = adata_pred[:, common_regions].X
    true_subset = adata_true[:, common_regions].X
    
    # Convert to dense if sparse
    pred_subset = np.asarray(pred_subset)
    true_subset = np.asarray(true_subset)
    
    # Rank transform for Spearman
    if method.lower() == "spearman":
        from scipy.stats import rankdata
        pred_subset = np.apply_along_axis(rankdata, 1, pred_subset)
        true_subset = np.apply_along_axis(rankdata, 1, true_subset)
    
    # Center and normalize each row
    def zscore_rows(X):
        X_mean = X.mean(axis=1, keepdims=True)
        X_std = X.std(axis=1, keepdims=True)
        X_std[X_std == 0] = 1  # Avoid division by zero
        return (X - X_mean) / X_std
    
    pred_z = zscore_rows(pred_subset)
    true_z = zscore_rows(true_subset)
    
    # Correlation matrix via matrix multiply
    # Each pred celltype (row of pred_z) vs each true celltype (row of true_z)
    n_regions = pred_z.shape[1]
    corr_matrix = (pred_z @ true_z.T) / (n_regions - 1)
    
    # Create DataFrame
    corr_df = pd.DataFrame(
        corr_matrix,
        index=adata_pred.obs_names,
        columns=adata_true.obs_names
    )
    
    return corr_df


def find_best_matches(
    corr_df: pd.DataFrame,
    top_k: int = 3
) -> pd.DataFrame:
    """
    For each predicted cell type, find best matching true cell types.
    
    Parameters
    ----------
    corr_df : DataFrame
        Correlation matrix from compare_across_celltypes()
    top_k : int
        Number of top matches to return
        
    Returns
    -------
    DataFrame
        Best matches for each predicted cell type
    """
    results = []
    
    for pred_ct in corr_df.index:
        correlations = corr_df.loc[pred_ct].sort_values(ascending=False)
        for rank, (true_ct, corr_val) in enumerate(correlations.head(top_k).items(), 1):
            results.append({
                'predicted_celltype': pred_ct,
                'true_celltype': true_ct,
                'correlation': corr_val,
                'rank': rank
            })
    
    return pd.DataFrame(results)


# ============================================================================
# COMPLETE WORKFLOW EXAMPLE
# ============================================================================

def minimal_workflow_example():
    """
    Complete minimal example: predict → compare → analyze.
    """
    
    # Load model (trained on species A, e.g., human)
    model = crested.load_model("path/to/human_model.h5")
    
    # Load data (from species B, e.g., chimp - DIFFERENT cell types)
    adata_chimp = ad.read_h5ad("path/to/chimp_data.h5ad")
    
    # Register genome if needed
    genome = crested.Genome("path/to/chimp_genome.fa")
    crested.register_genome(genome)
    
    print("="*60)
    print("STEP 1: Predict accessibility")
    print("="*60)
    
    # Predict - model will predict for its own cell types
    adata_pred = predict_accessibility(
        model=model,
        adata=adata_chimp,
        batch_size=16
    )
    
    print(f"\nPredictions shape: {adata_pred.shape}")
    print(f"Predicted cell types: {list(adata_pred.obs_names)}")
    print(f"Original cell types: {adata_pred.uns['input_cell_types']}")
    
    print("\n" + "="*60)
    print("STEP 2: Compare predicted vs true accessibility")
    print("="*60)
    
    # Compare across all cell type pairs
    corr_matrix = compare_across_celltypes(
        adata_pred=adata_pred,
        adata_true=adata_chimp,
        method="pearson"
    )
    
    print(f"\nCorrelation matrix shape: {corr_matrix.shape}")
    print("\nCorrelation matrix:")
    print(corr_matrix.round(3))
    
    print("\n" + "="*60)
    print("STEP 3: Find best matches")
    print("="*60)
    
    # Find which true cell types match which predicted cell types
    best_matches = find_best_matches(corr_matrix, top_k=3)
    print("\nBest matches:")
    print(best_matches.to_string(index=False))
    
    print("\n" + "="*60)
    print("STEP 4: Visualize")
    print("="*60)
    
    # Plot correlation matrix
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    plt.figure(figsize=(10, 8))
    sns.heatmap(
        corr_matrix,
        cmap='RdBu_r',
        center=0,
        vmin=-1,
        vmax=1,
        annot=True,
        fmt='.2f',
        square=True,
        cbar_kws={'label': 'Correlation'}
    )
    plt.title('Predicted Cell Types (rows) vs True Cell Types (cols)')
    plt.xlabel('True Cell Types (Chimp)')
    plt.ylabel('Predicted Cell Types (Human Model)')
    plt.tight_layout()
    plt.savefig('celltype_matching.pdf')
    print("Saved: celltype_matching.pdf")
    
    # Plot top correlations
    plt.figure(figsize=(12, 6))
    best_matches_top = best_matches[best_matches['rank'] == 1]
    plt.barh(
        best_matches_top['predicted_celltype'],
        best_matches_top['correlation']
    )
    plt.xlabel('Correlation')
    plt.title('Best Match for Each Predicted Cell Type')
    plt.tight_layout()
    plt.savefig('best_matches.pdf')
    print("Saved: best_matches.pdf")
    
    return adata_pred, corr_matrix, best_matches


# ============================================================================
# ULTRA-MINIMAL 10-LINE VERSION
# ============================================================================

def ultra_minimal_example():
    """
    Absolute minimal version - just 10 lines!
    """
    import crested
    import anndata as ad
    
    # Load
    model = crested.load_model("model.h5")
    adata = ad.read_h5ad("data.h5ad")
    
    # Predict
    adata_pred = predict_accessibility(model, adata)
    
    # Compare
    corr = compare_across_celltypes(adata_pred, adata)
    
    # Analyze
    print("Best matches:")
    print(find_best_matches(corr))


if __name__ == "__main__":
    print(__doc__)
    print("\nTo use:")
    print("1. Import: from minimal_predict import predict_accessibility, compare_across_celltypes")
    print("2. Predict: adata_pred = predict_accessibility(model, adata)")
    print("3. Compare: corr = compare_across_celltypes(adata_pred, adata)")
    print("4. Analyze: matches = find_best_matches(corr)")
    print("\nSee ultra_minimal_example() for complete 10-line workflow")
