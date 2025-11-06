"""
Example usage of crested_utils for cross-species predictions and comparisons.

This script demonstrates how to:
1. Load data for different species
2. Make predictions using a trained model
3. Align cell types across species
4. Compare predicted vs true values
5. Visualize correlations
"""

import anndata as ad
import crested
import matplotlib.pyplot as plt
import seaborn as sns
from crested_utils import (
    predict_regions,
    align_adata_cell_types,
    compare_predictions,
    rowwise_correlation,
    resize_region
)

# =============================================================================
# Example 1: Basic prediction on an AnnData object
# =============================================================================

def example_basic_prediction():
    """Make predictions on an existing AnnData object."""
    
    # Load your data
    adata = ad.read_h5ad("path/to/your/adata.h5ad")
    
    # Load your trained model
    model = crested.load_model("path/to/your/model")
    
    # Make predictions (automatically handles chunking for large datasets)
    adata_with_pred = predict_regions(
        model=model,
        regions=adata,  # Pass your AnnData directly
        batch_size=16,
        chunk_size=1000,  # Adjust based on your memory
        layer_name="predicted",
        verbose=True
    )
    
    # Now adata_with_pred has predictions in .layers['predicted']
    print(f"Shape: {adata_with_pred.shape}")
    print(f"Layers: {list(adata_with_pred.layers.keys())}")
    
    return adata_with_pred


# =============================================================================
# Example 2: Cross-species prediction with cell type alignment
# =============================================================================

def example_cross_species():
    """
    Make predictions for a different species and align cell types to a
    reference (e.g., human).
    """
    
    # Define reference cell types (e.g., from human)
    human_celltypes = [
        "T_cell", "B_cell", "NK_cell", "Macrophage", 
        "Dendritic_cell", "Epithelial", "Endothelial"
    ]
    
    # Load species data (e.g., chimpanzee)
    adata_chimp = ad.read_h5ad("chimpanzee_data.h5ad")
    
    # Load model (trained on human data)
    model = crested.load_model("human_model")
    
    # Register species genome
    genome_chimp = crested.Genome("path/to/panTro3/genome.fa")
    crested.register_genome(genome_chimp)
    
    # Make predictions WITH cell type alignment
    adata_chimp_pred = predict_regions(
        model=model,
        regions=adata_chimp,
        target_cell_types=human_celltypes,
        align_cell_types=True,  # Align to human cell types
        fill_missing_with_zeros=True,  # Fill missing cell types with zeros
        chunk_size=2000,
        verbose=True
    )
    
    # Now all species have the same cell types in the same order
    print(f"Aligned cell types: {list(adata_chimp_pred.obs_names)}")
    print(f"All match human: {list(adata_chimp_pred.obs_names) == human_celltypes}")
    
    return adata_chimp_pred


# =============================================================================
# Example 3: Compare predictions vs true values
# =============================================================================

def example_compare_predictions():
    """Compare predicted vs true values and visualize."""
    
    # Load data with predictions (from previous examples)
    adata = ad.read_h5ad("adata_with_predictions.h5ad")
    
    # Compare predictions
    corr_matrix, self_corr = compare_predictions(
        adata,
        prediction_layer="predicted",
        true_layer=None,  # None means use .X
        method="pearson"
    )
    
    print("\nSelf-correlation per cell type:")
    print(self_corr.sort_values('correlation', ascending=False))
    print(f"\nMean correlation: {self_corr['correlation'].mean():.3f}")
    
    # Visualize correlation matrix
    plt.figure(figsize=(12, 10))
    sns.heatmap(
        corr_matrix,
        cmap='RdBu_r',
        center=0,
        vmin=-1,
        vmax=1,
        square=True,
        cbar_kws={'label': 'Correlation'}
    )
    plt.title('Predicted vs True: Row-wise Correlation')
    plt.xlabel('Predicted')
    plt.ylabel('True')
    plt.tight_layout()
    plt.savefig('correlation_heatmap.pdf')
    
    return corr_matrix, self_corr


# =============================================================================
# Example 4: Full cross-species analysis pipeline
# =============================================================================

def example_full_pipeline():
    """
    Complete pipeline for cross-species prediction and comparison.
    """
    
    # Configuration
    species_list = ["chimpanzee", "gorilla", "macaque"]
    genome_map = {
        "chimpanzee": "panTro3",
        "gorilla": "gorGor4",
        "macaque": "Mmul10"
    }
    base_genome_path = "/path/to/genomes/reference_"
    
    # Load human model and define reference cell types
    model = crested.load_model("human_model")
    human_celltypes = [
        "T_cell", "B_cell", "NK_cell", "Macrophage", 
        "Dendritic_cell", "Epithelial", "Endothelial"
    ]
    
    # Store results
    results = {}
    
    for species in species_list:
        print(f"\n{'='*60}")
        print(f"Processing {species}")
        print(f"{'='*60}")
        
        # Load species data
        adata_species = ad.read_h5ad(f"{species}_data.h5ad")
        
        # Register genome
        genome_name = genome_map[species]
        genome_path = f"{base_genome_path}/{genome_name}/fasta/genome.fa"
        genome = crested.Genome(genome_path)
        crested.register_genome(genome)
        
        # Optional: Filter to most specific regions for this species
        # This can speed up computation if you have many regions
        common_celltypes = [ct for ct in human_celltypes if ct in adata_species.obs_names]
        adata_filtered = adata_species[common_celltypes, :].copy()
        
        # You can also filter to top-k most specific regions
        # crested.pp.sort_and_filter_regions_on_specificity(
        #     adata_filtered, top_k=500, method="proportion"
        # )
        
        # Make predictions with alignment
        adata_pred = predict_regions(
            model=model,
            regions=adata_filtered,
            target_cell_types=human_celltypes,
            align_cell_types=True,
            fill_missing_with_zeros=True,
            chunk_size=2000,
            disk_based_saving=False,  # Use True for very large datasets
            verbose=True
        )
        
        # Compare predictions for common cell types only
        corr_matrix, self_corr = compare_predictions(
            adata_pred,
            prediction_layer="predicted",
            cell_types=common_celltypes,  # Only compare cell types present in species
            method="pearson"
        )
        
        # Store results
        results[species] = {
            'adata': adata_pred,
            'correlation_matrix': corr_matrix,
            'self_correlation': self_corr,
            'mean_correlation': self_corr['correlation'].mean()
        }
        
        print(f"\n{species} results:")
        print(f"  Mean correlation: {results[species]['mean_correlation']:.3f}")
        print(f"  Common cell types: {len(common_celltypes)}")
    
    # Visualize comparison across species
    fig, axes = plt.subplots(1, len(species_list), figsize=(6*len(species_list), 5))
    
    for idx, species in enumerate(species_list):
        ax = axes[idx] if len(species_list) > 1 else axes
        corr_mat = results[species]['correlation_matrix']
        
        sns.heatmap(
            corr_mat,
            ax=ax,
            cmap='RdBu_r',
            center=0,
            vmin=-1,
            vmax=1,
            square=True,
            cbar_kws={'label': 'Correlation'}
        )
        
        mean_corr = results[species]['mean_correlation']
        ax.set_title(f'{species.capitalize()}\nMean corr: {mean_corr:.3f}')
        ax.set_xlabel('Predicted')
        ax.set_ylabel('True')
    
    plt.tight_layout()
    plt.savefig('cross_species_comparison.pdf')
    
    return results


# =============================================================================
# Example 5: Manual cell type alignment (if needed separately)
# =============================================================================

def example_manual_alignment():
    """Manually align cell types if you already have predictions."""
    
    # Load data with predictions but mismatched cell types
    adata_species = ad.read_h5ad("species_with_predictions.h5ad")
    
    # Define target cell types
    target_celltypes = [
        "T_cell", "B_cell", "NK_cell", "Macrophage", 
        "Dendritic_cell", "Epithelial", "Endothelial"
    ]
    
    # Align (this preserves layers including predictions)
    adata_aligned = align_adata_cell_types(
        adata_species,
        target_cell_types=target_celltypes,
        fill_missing=True,  # True = add zeros for missing, False = drop missing
        verbose=True
    )
    
    print(f"Original shape: {adata_species.shape}")
    print(f"Aligned shape: {adata_aligned.shape}")
    print(f"Cell types match target: {list(adata_aligned.obs_names) == target_celltypes}")
    
    return adata_aligned


# =============================================================================
# Example 6: Region resizing (for contribution scores, etc.)
# =============================================================================

def example_region_resizing():
    """Resize regions for analysis like contribution scores."""
    
    # Original regions
    regions_of_interest = [
        "chr21:18402602-18404716",  # Example: FIRE enhancer
        "chr12:112924652-112925561"
    ]
    
    # Resize to fixed length (e.g., for model input)
    target_length = 2114
    regions_resized = [resize_region(r, target_length) for r in regions_of_interest]
    
    print("Original regions:")
    for r in regions_of_interest:
        print(f"  {r}")
    
    print(f"\nResized to {target_length}bp:")
    for r in regions_resized:
        print(f"  {r}")
    
    # Use resized regions for contribution scores
    # model = crested.load_model("your_model")
    # adata = ad.read_h5ad("your_data.h5ad")
    # class_idx = list(adata.obs_names.get_indexer(classes_of_interest))
    # scores, sequences = crested.tl.contribution_scores(
    #     regions_resized,
    #     target_idx=class_idx,
    #     model=model,
    #     batch_size=24
    # )
    
    return regions_resized


# =============================================================================
# Example 7: Direct row-wise correlation (advanced)
# =============================================================================

def example_direct_correlation():
    """
    Compute correlation directly between two DataFrames.
    Useful for custom analyses.
    """
    import pandas as pd
    
    # Load predictions from different models/species
    adata1 = ad.read_h5ad("model1_predictions.h5ad")
    adata2 = ad.read_h5ad("model2_predictions.h5ad")
    
    # Convert to DataFrames
    pred1 = pd.DataFrame(
        adata1.layers['predicted'],
        index=adata1.obs_names,
        columns=adata1.var_names
    )
    pred2 = pd.DataFrame(
        adata2.layers['predicted'],
        index=adata2.obs_names,
        columns=adata2.var_names
    )
    
    # Compute correlation
    corr = rowwise_correlation(pred1, pred2, method="pearson")
    
    print(f"Correlation matrix shape: {corr.shape}")
    print(f"Diagonal (self-correlation):")
    print(pd.Series(corr.values.diagonal(), index=corr.index))
    
    return corr


# =============================================================================
# Main execution
# =============================================================================

if __name__ == "__main__":
    
    print("CREsted Utils Examples")
    print("=" * 60)
    print("\nThese examples demonstrate the main functions.")
    print("Uncomment the example you want to run.\n")
    
    # Example 1: Basic prediction
    # adata_pred = example_basic_prediction()
    
    # Example 2: Cross-species with alignment
    # adata_species_pred = example_cross_species()
    
    # Example 3: Compare predictions
    # corr_mat, self_corr = example_compare_predictions()
    
    # Example 4: Full pipeline
    # results = example_full_pipeline()
    
    # Example 5: Manual alignment
    # adata_aligned = example_manual_alignment()
    
    # Example 6: Region resizing
    # regions_resized = example_region_resizing()
    
    # Example 7: Direct correlation
    # corr = example_direct_correlation()
    
    print("\nDone! Uncomment examples in the main block to run them.")
