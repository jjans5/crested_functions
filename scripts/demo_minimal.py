"""
Quick demo of minimal_predict.py with synthetic data.

This shows the workflow without needing real data files.
"""

import numpy as np
import pandas as pd
import anndata as ad
from minimal_predict import (
    predict_accessibility,
    compare_across_celltypes,
    find_best_matches
)


def create_synthetic_model():
    """Create a fake model for demo purposes."""
    class FakeModel:
        """Fake model that returns random predictions."""
        def __init__(self, n_celltypes=3):
            self.n_celltypes = n_celltypes
    
    return FakeModel(n_celltypes=3)


def create_synthetic_adata(n_celltypes=5, n_regions=100, cell_names=None):
    """Create synthetic AnnData for demo."""
    np.random.seed(42)
    
    # Create data
    X = np.random.rand(n_celltypes, n_regions).astype(np.float32)
    
    # Cell type names
    if cell_names is None:
        cell_names = [f'CellType_{i}' for i in range(n_celltypes)]
    
    # Region names
    region_names = [f'chr1:{i*1000}-{i*1000+500}' for i in range(n_regions)]
    
    # Create AnnData
    obs = pd.DataFrame(index=cell_names)
    var = pd.DataFrame(index=region_names)
    
    adata = ad.AnnData(X=X, obs=obs, var=var)
    
    return adata


# Monkey-patch crested.tl.predict for demo
def fake_predict(adata, model, batch_size=16):
    """Fake prediction function."""
    n_regions = adata.n_vars
    n_celltypes = model.n_celltypes
    
    # Return random predictions
    np.random.seed(123)
    predictions = np.random.rand(n_regions, n_celltypes).astype(np.float32)
    
    return predictions


def demo_minimal_workflow():
    """
    Demonstrate the complete workflow with synthetic data.
    
    Scenario: 
    - Model trained on human: T_cell, B_cell, NK_cell (3 types)
    - Data from chimp: CD4_T, CD8_T, B_cell, Monocyte, Dendritic (5 types)
    - ZERO overlap in cell type names!
    """
    
    print("="*70)
    print("DEMO: Predict & Compare with Zero Cell Type Overlap")
    print("="*70)
    
    # Create synthetic data
    print("\nCreating synthetic data...")
    
    # "Human" model (3 cell types)
    model = create_synthetic_model()
    print(f"✓ Model: {model.n_celltypes} cell types (human T_cell, B_cell, NK_cell)")
    
    # "Chimp" data (5 different cell types)
    chimp_celltypes = ['CD4_T', 'CD8_T', 'B_cell_chimp', 'Monocyte', 'Dendritic']
    adata_chimp = create_synthetic_adata(
        n_celltypes=5,
        n_regions=100,
        cell_names=chimp_celltypes
    )
    print(f"✓ Chimp data: {adata_chimp.n_obs} cell types × {adata_chimp.n_vars} regions")
    print(f"  Cell types: {list(adata_chimp.obs_names)}")
    
    # Monkey-patch crested for demo
    import crested
    if not hasattr(crested, 'tl'):
        class TL:
            predict = staticmethod(fake_predict)
        crested.tl = TL()
    else:
        crested.tl.predict = fake_predict
    
    print("\n" + "="*70)
    print("STEP 1: Predict accessibility")
    print("="*70)
    
    # Predict - model predicts for its own cell types
    adata_pred = predict_accessibility(model, adata_chimp)
    
    print(f"\n✓ Predictions created!")
    print(f"  Shape: {adata_pred.shape}")
    print(f"  Model's cell types: {list(adata_pred.obs_names)}")
    print(f"  Regions: {adata_pred.n_vars}")
    
    print("\n" + "="*70)
    print("STEP 2: Compare predicted vs true")
    print("="*70)
    
    # Compare all predicted vs all true cell types
    corr_matrix = compare_across_celltypes(
        adata_pred=adata_pred,
        adata_true=adata_chimp,
        method="pearson"
    )
    
    print(f"\n✓ Correlation matrix computed!")
    print(f"  Shape: {corr_matrix.shape}")
    print(f"  (rows = model's {adata_pred.n_obs} celltypes, cols = true {adata_chimp.n_obs} celltypes)")
    print("\nCorrelation Matrix:")
    print(corr_matrix.round(3).to_string())
    
    print("\n" + "="*70)
    print("STEP 3: Find best matches")
    print("="*70)
    
    # Find which true cell types match predicted ones
    best_matches = find_best_matches(corr_matrix, top_k=3)
    
    print("\n✓ Best matches identified!")
    print("\nTop 3 matches for each predicted cell type:")
    print(best_matches.to_string(index=False))
    
    # Summary
    print("\n" + "="*70)
    print("INTERPRETATION")
    print("="*70)
    print("\nThe correlation matrix shows:")
    print("- Each row = one of the model's predicted cell types")
    print("- Each column = one of the true cell types in the data")
    print("- High correlation = similar accessibility patterns")
    print("\nEven though cell type NAMES don't match,")
    print("you can identify which cell types are functionally similar!")
    
    print("\n" + "="*70)
    print("YOUR 3-LINE MINIMAL CODE")
    print("="*70)
    print("""
    from minimal_predict import predict_accessibility, compare_across_celltypes, find_best_matches
    
    adata_pred = predict_accessibility(model, adata)
    corr = compare_across_celltypes(adata_pred, adata)
    matches = find_best_matches(corr)
    """)
    
    return adata_pred, corr_matrix, best_matches


if __name__ == "__main__":
    # Run demo
    adata_pred, corr, matches = demo_minimal_workflow()
    
    print("\n" + "="*70)
    print("DEMO COMPLETE!")
    print("="*70)
    print("\nNow try with your real data:")
    print("  1. Load your model: model = crested.load_model('model.h5')")
    print("  2. Load your data: adata = ad.read_h5ad('data.h5ad')")
    print("  3. Run the 3 lines above!")
