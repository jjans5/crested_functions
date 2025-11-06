"""
Simple tests to verify crested_utils functions work correctly.

Run this after installing dependencies to make sure everything is working.
"""

import numpy as np
import pandas as pd
import anndata as ad
from crested_utils import (
    align_adata_cell_types,
    rowwise_correlation,
    resize_region,
    compare_predictions
)


def test_resize_region():
    """Test region resizing."""
    print("Testing resize_region()...")
    
    # Test basic resizing
    result = resize_region("chr1:1000-2000", 500)
    assert "chr1:" in result
    start, end = map(int, result.split(":")[1].split("-"))
    assert end - start == 500
    print("  ✓ Basic resizing works")
    
    # Test with summit
    result = resize_region("chr1:1000-2000", 500, summit=1600)
    start, end = map(int, result.split(":")[1].split("-"))
    assert end - start == 500
    assert start == 1350  # 1600 - 250
    print("  ✓ Summit-based resizing works")
    
    print("✅ resize_region() passed all tests\n")


def test_rowwise_correlation():
    """Test row-wise correlation."""
    print("Testing rowwise_correlation()...")
    
    # Create test data
    np.random.seed(42)
    data1 = np.random.rand(5, 100)
    data2 = np.random.rand(5, 100)
    
    df1 = pd.DataFrame(data1, index=['A', 'B', 'C', 'D', 'E'])
    df2 = pd.DataFrame(data2, index=['A', 'B', 'C', 'D', 'E'])
    
    # Test Pearson
    corr = rowwise_correlation(df1, df2, method='pearson')
    assert corr.shape == (5, 5)
    assert -1 <= corr.values.min() <= 1
    assert -1 <= corr.values.max() <= 1
    print("  ✓ Pearson correlation works")
    
    # Test Spearman
    corr = rowwise_correlation(df1, df2, method='spearman')
    assert corr.shape == (5, 5)
    print("  ✓ Spearman correlation works")
    
    # Test self-correlation
    corr = rowwise_correlation(df1, df1, method='pearson')
    diag = np.diag(corr.values)
    assert np.allclose(diag, 1.0)  # Self-correlation should be 1
    print("  ✓ Self-correlation is 1.0")
    
    print("✅ rowwise_correlation() passed all tests\n")


def test_align_adata_cell_types():
    """Test cell type alignment."""
    print("Testing align_adata_cell_types()...")
    
    # Create test AnnData
    n_cells = 5
    n_regions = 100
    
    X = np.random.rand(n_cells, n_regions).astype(np.float32)
    obs = pd.DataFrame(index=['T_cell', 'B_cell', 'NK_cell', 'Macrophage', 'Dendritic'])
    var = pd.DataFrame(index=[f'region_{i}' for i in range(n_regions)])
    
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.layers['predictions'] = X * 1.1
    adata.obsm['weights'] = np.random.rand(n_cells, 10)
    
    # Test with all cell types present
    target = ['T_cell', 'B_cell', 'NK_cell', 'Macrophage', 'Dendritic']
    adata_aligned = align_adata_cell_types(adata, target, verbose=False)
    assert list(adata_aligned.obs_names) == target
    assert adata_aligned.shape == adata.shape
    print("  ✓ Alignment with all cell types present works")
    
    # Test with missing cell types (fill with zeros)
    target = ['T_cell', 'B_cell', 'NK_cell', 'Macrophage', 'Dendritic', 'Epithelial', 'Endothelial']
    adata_aligned = align_adata_cell_types(adata, target, fill_missing=True, verbose=False)
    assert list(adata_aligned.obs_names) == target
    assert adata_aligned.shape == (7, n_regions)
    # Check that added cell types are zeros
    assert np.allclose(adata_aligned.X[5, :], 0)  # Epithelial
    assert np.allclose(adata_aligned.X[6, :], 0)  # Endothelial
    print("  ✓ Alignment with missing cell types (fill zeros) works")
    
    # Test with missing cell types (drop them)
    target = ['T_cell', 'B_cell', 'Macrophage', 'Epithelial']
    adata_aligned = align_adata_cell_types(adata, target, fill_missing=False, verbose=False)
    assert list(adata_aligned.obs_names) == ['T_cell', 'B_cell', 'Macrophage']
    assert adata_aligned.shape == (3, n_regions)
    print("  ✓ Alignment with missing cell types (drop) works")
    
    # Test layer preservation
    target = ['T_cell', 'B_cell', 'NK_cell', 'Macrophage', 'Dendritic']
    adata_aligned = align_adata_cell_types(adata, target, verbose=False)
    assert 'predictions' in adata_aligned.layers
    assert 'weights' in adata_aligned.obsm
    print("  ✓ Layers and obsm preserved")
    
    # Test dtype preservation
    assert adata_aligned.X.dtype == np.float32
    print("  ✓ dtype preserved (float32)")
    
    print("✅ align_adata_cell_types() passed all tests\n")


def test_compare_predictions():
    """Test prediction comparison."""
    print("Testing compare_predictions()...")
    
    # Create test AnnData with predictions
    n_cells = 5
    n_regions = 100
    
    X_true = np.random.rand(n_cells, n_regions).astype(np.float32)
    # Predictions correlated with true + noise
    X_pred = X_true * 0.8 + np.random.rand(n_cells, n_regions).astype(np.float32) * 0.2
    
    obs = pd.DataFrame(index=['T_cell', 'B_cell', 'NK_cell', 'Macrophage', 'Dendritic'])
    var = pd.DataFrame(index=[f'region_{i}' for i in range(n_regions)])
    
    adata = ad.AnnData(X=X_true, obs=obs, var=var)
    adata.layers['predicted'] = X_pred
    
    # Test comparison
    corr_matrix, self_corr = compare_predictions(adata, prediction_layer='predicted')
    
    assert corr_matrix.shape == (n_cells, n_cells)
    assert len(self_corr) == n_cells
    assert 'cell_type' in self_corr.columns
    assert 'correlation' in self_corr.columns
    print("  ✓ Comparison returns correct shapes")
    
    # Check correlations are reasonable
    assert self_corr['correlation'].min() > 0  # Should be positive since pred ~ true
    assert self_corr['correlation'].max() <= 1.0
    print("  ✓ Correlations are in reasonable range")
    
    # Test with subset of cell types
    corr_matrix, self_corr = compare_predictions(
        adata,
        cell_types=['T_cell', 'B_cell'],
        prediction_layer='predicted'
    )
    assert corr_matrix.shape == (2, 2)
    assert len(self_corr) == 2
    print("  ✓ Subsetting cell types works")
    
    # Test Spearman
    corr_matrix, self_corr = compare_predictions(
        adata,
        method='spearman'
    )
    assert corr_matrix.shape == (n_cells, n_cells)
    print("  ✓ Spearman correlation works")
    
    print("✅ compare_predictions() passed all tests\n")


def run_all_tests():
    """Run all tests."""
    print("="*60)
    print("Running crested_utils tests")
    print("="*60 + "\n")
    
    try:
        test_resize_region()
        test_rowwise_correlation()
        test_align_adata_cell_types()
        test_compare_predictions()
        
        print("="*60)
        print("🎉 All tests passed!")
        print("="*60)
        print("\nYour crested_utils module is working correctly.")
        print("You can now use it for your analysis.")
        
    except AssertionError as e:
        print(f"\n❌ Test failed: {e}")
        raise
    except Exception as e:
        print(f"\n❌ Unexpected error: {e}")
        raise


if __name__ == "__main__":
    run_all_tests()
