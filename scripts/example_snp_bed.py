"""
Example usage of SNP mutagenesis from BED files.

This script demonstrates how to analyze SNPs from a BED file.
"""

import numpy as np
import pandas as pd
import anndata as ad
from src.insilico_mutagenesis_vect import snp_mutagenesis_from_bed
import crested


def create_example_bed_file(filepath="example_snps.bed"):
    """Create an example BED file with SNPs."""
    # Example 1: BED file with positions only (will test all mutations)
    bed_positions_only = """chr1\t1000\t1001
chr1\t2000\t2001
chr2\t3000\t3001"""
    
    with open(filepath, 'w') as f:
        f.write(bed_positions_only)
    
    print(f"Created example BED file: {filepath}")
    print("Format: chrom, start, end (will test all mutations)")
    
    # Example 2: BED file with ref/alt alleles
    bed_with_alleles_path = filepath.replace('.bed', '_with_alleles.bed')
    bed_with_alleles = """chr1\t1000\t1001\tA\tG
chr1\t2000\t2001\tC\tT
chr2\t3000\t3001\tG\tA"""
    
    with open(bed_with_alleles_path, 'w') as f:
        f.write(bed_with_alleles)
    
    print(f"Created example BED file with alleles: {bed_with_alleles_path}")
    print("Format: chrom, start, end, ref, alt")


def example_snp_analysis_from_bed():
    """Example: Analyze SNPs from a BED file."""
    
    print("\n" + "="*60)
    print("Example: SNP Analysis from BED File")
    print("="*60 + "\n")
    
    # Load your model and data
    model_path = "path/to/your/model.keras"
    genome_path = "path/to/genome"
    
    # Load model
    print("Loading model...")
    model = crested.load_model(model_path)
    
    # Load genome
    print("Loading genome...")
    genome = crested.Genome(genome_path)
    
    # Create example AnnData with cell type names
    cell_types = ["CD4_T", "CD8_T", "B_cell", "Monocyte"]
    adata = ad.AnnData(
        X=np.zeros((len(cell_types), 1)),
        obs=pd.DataFrame(index=cell_types)
    )
    
    # Create example BED files
    create_example_bed_file("example_snps.bed")
    
    print("\n" + "-"*60)
    print("Example 1: BED file with positions only")
    print("-"*60)
    
    # Analyze SNPs (will test all possible mutations)
    results = snp_mutagenesis_from_bed(
        bed_file="example_snps.bed",
        model=model,
        adata=adata,
        genome=genome,
        seq_length=2114,  # Model input length
        batch_size=8
    )
    
    print(f"\nResults shape: {results.shape}")
    print(f"Columns: {list(results.columns)}")
    print("\nFirst few results:")
    print(results.head())
    
    print("\n" + "-"*60)
    print("Example 2: BED file with ref/alt specified")
    print("-"*60)
    
    # Analyze SNPs with specific ref/alt
    results_specific = snp_mutagenesis_from_bed(
        bed_file="example_snps_with_alleles.bed",
        model=model,
        adata=adata,
        genome=genome
    )
    
    print(f"\nResults shape: {results_specific.shape}")
    print("\nFirst few results:")
    print(results_specific.head())
    
    # Analyze results
    print("\n" + "-"*60)
    print("Analysis Summary")
    print("-"*60)
    
    # Find SNPs with largest effects
    print("\nTop 10 SNPs by absolute log2fc:")
    top_snps = results.nlargest(10, 'log2fc')[
        ['mut_id', 'cell_type', 'log2fc', 'delta']
    ]
    print(top_snps)
    
    # Summary by cell type
    print("\nMean absolute log2fc by cell type:")
    cell_type_summary = results.groupby('cell_type')['log2fc'].apply(
        lambda x: np.mean(np.abs(x))
    ).sort_values(ascending=False)
    print(cell_type_summary)
    
    # Save results
    output_file = "snp_analysis_results.tsv"
    results.to_csv(output_file, sep='\t', index=False)
    print(f"\n✅ Results saved to: {output_file}")


def example_filter_significant_snps():
    """Example: Filter SNPs by effect size."""
    
    print("\n" + "="*60)
    print("Example: Filter Significant SNPs")
    print("="*60 + "\n")
    
    # Assuming you have results from previous analysis
    # results = pd.read_csv("snp_analysis_results.tsv", sep='\t')
    
    # Filter SNPs with log2fc > 1 (2-fold change)
    # significant = results[np.abs(results['log2fc']) > 1]
    
    # Group by SNP to see effect across cell types
    # snp_effects = significant.groupby('mut_id').agg({
    #     'log2fc': ['mean', 'std', 'min', 'max'],
    #     'cell_type': 'count'
    # })
    
    print("Example filtering code:")
    print("""
    # Load results
    results = pd.read_csv("snp_analysis_results.tsv", sep='\\t')
    
    # Filter by effect size (log2fc > 1 means 2-fold change)
    significant = results[np.abs(results['log2fc']) > 1]
    
    # Find SNPs with consistent effects across cell types
    snp_summary = significant.groupby('mut_id').agg({
        'log2fc': ['mean', 'std', 'min', 'max'],
        'cell_type': 'count'
    })
    
    # SNPs affecting multiple cell types
    multi_cell_type = snp_summary[snp_summary[('cell_type', 'count')] > 2]
    print(multi_cell_type)
    """)


if __name__ == "__main__":
    print("SNP Mutagenesis from BED File - Example Usage")
    print("=" * 60)
    
    print("\nNote: This example requires:")
    print("- A trained CREsted model")
    print("- A reference genome")
    print("- A BED file with SNP positions")
    
    print("\nTo run this example, update the paths in the code:")
    print("- model_path: path to your .keras model file")
    print("- genome_path: path to your genome directory")
    print("- bed_file: path to your BED file with SNPs")
    
    # Uncomment to run examples (after updating paths):
    # example_snp_analysis_from_bed()
    # example_filter_significant_snps()
    
    print("\n✅ See function examples above for usage patterns")
