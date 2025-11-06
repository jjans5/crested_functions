import numpy as np
import pandas as pd
from typing import Optional, Union, List


def insilico_mutagenesis_vect(
    seq: str,
    model,
    cell_types: Optional[Union[List[str], 'anndata.AnnData']] = None,
    chrom: str | None = None,
    start: int | None = None,   # genomic start for seq[0]; if None we report 0-based pos0
    skip_ambiguous: bool = True,
    batch_size=2,
    return_long: bool = True,   # also return tidy long table with deltas
):
    """
    Vectorized single-nucleotide in-silico mutagenesis (ISM).

    Parameters
    ----------
    seq : str
        DNA sequence to mutate
    model : CREsted model
        Trained model for predictions
    cell_types : list of str or AnnData, optional
        Cell type names for predictions. Can be:
        - List of strings: ['CellType1', 'CellType2', ...]
        - AnnData object: will use adata.obs_names
        - None: will auto-generate names from model output
    chrom : str, optional
        Chromosome name for annotation
    start : int, optional
        Genomic start position for seq[0]
    skip_ambiguous : bool, default=True
        Skip positions with ambiguous nucleotides
    batch_size : int, default=2
        Batch size for predictions
    return_long : bool, default=True
        Return tidy long-format DataFrame with log2fc

    Returns
    -------
    dict with:
      wt_pred:      pd.Series (index=cell_types)
      mut_df_wide:  pd.DataFrame with metadata columns ['chrom','genomic_pos','pos0','ref','alt','mut_id']
                    plus one column per cell type (absolute predictions)
      mut_df_long:  (if return_long) tidy DataFrame with columns
                    ['chrom','genomic_pos','pos0','ref','alt','mut_id','cell_type','value','log2fc','delta']
                    where log2fc is log2(mutant/wildtype) and delta is (mutant - wildtype)
    """
    import crested

    seq = seq.upper()
    bases = ("A", "C", "G", "T")
    
    # Handle cell_types input
    if cell_types is None:
        # Will infer from model output shape
        cell_types_list = None
    elif hasattr(cell_types, 'obs_names'):
        # It's an AnnData object
        cell_types_list = list(map(str, cell_types.obs_names))
    elif isinstance(cell_types, (list, tuple)):
        # It's already a list
        cell_types_list = list(map(str, cell_types))
    else:
        raise TypeError("cell_types must be a list of strings, an AnnData object, or None")
    
    L = len(seq)

    # 1) WT prediction
    wt = crested.tl.predict(input=[seq], model=model, batch_size=batch_size)
    wt = np.asarray(wt, dtype=np.float32).reshape(1, -1)
    
    # Infer cell types from output if not provided
    if cell_types_list is None:
        C = wt.shape[1]
        cell_types_list = [f"CellType_{i}" for i in range(C)]
    else:
        C = len(cell_types_list)
        if wt.shape[1] != C:
            raise ValueError(
                f"Prediction size mismatch: got {wt.shape[1]} cell types from model, "
                f"expected {C} from cell_types parameter"
            )
    
    wt_pred = pd.Series(wt[0], index=cell_types_list, name="wildtype")

    # 2) Build all mutants (vectorized): 3 mutants per valid position
    mutant_seqs = []
    chrom_list, gpos_list, pos0_list, ref_list, alt_list = [], [], [], [], []

    for i, ref in enumerate(seq):
        if ref not in bases:
            if skip_ambiguous:
                continue
            else:
                # If not skipping, we still only allow A/C/G/T as alts (no N->N)
                pass
        for alt in bases:
            if alt == ref:
                continue
            # assemble mutated sequence
            mutant_seqs.append(seq[:i] + alt + seq[i+1:])
            chrom_list.append(chrom)
            gpos_list.append((start + i) if start is not None else None)
            pos0_list.append(i)
            ref_list.append(ref)
            alt_list.append(alt)

    if not mutant_seqs:
        mut_df_wide = pd.DataFrame(columns=["chrom","genomic_pos","pos0","ref","alt","mut_id"] + cell_types_list)
        out = {"wt_pred": wt_pred, "mut_df_wide": mut_df_wide}
        if return_long:
            out["mut_df_long"] = pd.DataFrame(columns=["chrom","genomic_pos","pos0","ref","alt","mut_id","cell_type","value","delta"])
        return out

    # 3) Predict all mutants in a single call
    preds = crested.tl.predict(input=mutant_seqs, model=model, batch_size=batch_size)
    preds = np.asarray(preds, dtype=np.float32)
    if preds.shape != (len(mutant_seqs), C):
        raise ValueError(f"Predictions shape {preds.shape} does not match (n_mutants, n_celltypes)=({len(mutant_seqs)}, {C})")

    # 4) Metadata + wide table
    meta_df = pd.DataFrame({
        "chrom": chrom_list,
        "genomic_pos": gpos_list,
        "pos0": pos0_list,
        "ref": ref_list,
        "alt": alt_list,
    })
    # Stable, human-readable mutation identifier
    if start is not None and chrom is not None:
        meta_df["mut_id"] = (
            meta_df["chrom"].astype(str) + ":" +
            meta_df["genomic_pos"].astype(int).astype(str) + ":" +
            meta_df["ref"] + ">" + meta_df["alt"]
        )
    elif start is not None:
        meta_df["mut_id"] = (
            meta_df["genomic_pos"].astype(int).astype(str) + ":" +
            meta_df["ref"] + ">" + meta_df["alt"]
        )
    else:
        meta_df["mut_id"] = (
            meta_df["pos0"].astype(int).astype(str) + ":" +
            meta_df["ref"] + ">" + meta_df["alt"]
        )

    pred_df = pd.DataFrame(preds, columns=cell_types_list)
    mut_df_wide = pd.concat([meta_df, pred_df], axis=1)

    out = {"wt_pred": wt_pred, "mut_df_wide": mut_df_wide}

    # 5) Optional: tidy long table with deltas (log2 fold change)
    if return_long:
        long = mut_df_wide.melt(
            id_vars=["chrom","genomic_pos","pos0","ref","alt","mut_id"],
            value_vars=cell_types,
            var_name="cell_type",
            value_name="value"
        )
        # map wt once (fast)
        wt_map = wt_pred.to_dict()
        wt_values = long["cell_type"].map(wt_map)
        
        # Calculate log2 fold change: log2(mut/wt)
        # Add small epsilon to avoid division by zero
        epsilon = 1e-10
        long["log2fc"] = np.log2((long["value"] + epsilon) / (wt_values + epsilon)).astype(np.float32)
        
        # Also keep absolute delta for reference
        long["delta"] = (long["value"] - wt_values).astype(np.float32)
        
        out["mut_df_long"] = long

    return out


def snp_mutagenesis_from_bed(
    bed_file: Union[str, pd.DataFrame],
    model,
    genome,
    cell_types: Optional[Union[List[str], 'anndata.AnnData']] = None,
    seq_length: Optional[int] = None,
    batch_size: int = 32,  # Increased default for speed
    return_long: bool = True,
    skip_ambiguous: bool = True,
    chunk_size: int = 1000,  # Process SNPs in chunks to avoid memory issues
) -> pd.DataFrame:
    """
    Perform in-silico mutagenesis for SNPs from a BED file or DataFrame.
    
    This function:
    1. Reads SNP positions from a BED file or accepts a DataFrame
    2. Extracts sequences centered on each SNP (using model's seq_length)
    3. Predicts accessibility for wildtype and mutant sequences
    4. Returns log2 fold changes for each SNP-cell_type combination
    
    Parameters
    ----------
    bed_file : str or pd.DataFrame
        Path to BED file or a pandas DataFrame. Can have formats:
        - 3 columns: chrom, start, end (will test all mutations)
        - 4+ columns: chrom, start, end, name, ref, alt (will test specific ref>alt)
        - 5+ columns: chrom, start, end, ref, alt (name optional)
        - Flexible: ref and alt can be in any columns after position 3
        
        If DataFrame, should have columns matching the BED format above.
    model : CREsted model
        Trained model for predictions
    genome : crested.Genome
        Genome object for sequence extraction
    cell_types : list of str or AnnData, optional
        Cell type names for predictions. Can be:
        - List of strings: ['CellType1', 'CellType2', ...]
        - AnnData object: will use adata.obs_names
        - None: will auto-generate names from model output
    seq_length : int, optional
        Sequence length for model input. If None, inferred from model.
    batch_size : int, default=32
        Batch size for predictions (larger = faster but more memory)
    return_long : bool, default=True
        Return long-format DataFrame with log2fc values
    skip_ambiguous : bool, default=True
        Skip positions with ambiguous nucleotides (N, etc.)
    chunk_size : int, default=1000
        Number of SNPs to process per chunk (to avoid memory issues with large datasets)
        
    Returns
    -------
    pd.DataFrame
        Long-format DataFrame with columns:
        ['chrom', 'pos', 'ref', 'alt', 'mut_id', 'cell_type', 'wt_pred', 'mut_pred', 'log2fc', 'delta']
        
    Examples
    --------
    >>> # From BED file with just positions (will test all mutations)
    >>> results = snp_mutagenesis_from_bed(
    ...     bed_file="snps.bed",
    ...     model=model,
    ...     genome=genome,
    ...     cell_types=['CellType1', 'CellType2']  # or adata, or None
    ... )
    
    >>> # From BED file with ref/alt specified
    >>> results = snp_mutagenesis_from_bed(
    ...     bed_file="snps_with_alleles.bed",
    ...     model=model,
    ...     genome=genome,
    ...     cell_types=adata  # Can pass AnnData directly
    ... )
    
    >>> # From DataFrame with auto-generated cell type names
    >>> import pandas as pd
    >>> snps_df = pd.DataFrame({
    ...     'chrom': ['chr1', 'chr2'],
    ...     'start': [1000, 2000],
    ...     'end': [1001, 2001],
    ...     'ref': ['A', 'C'],
    ...     'alt': ['G', 'T']
    ... })
    >>> results = snp_mutagenesis_from_bed(
    ...     bed_file=snps_df,
    ...     model=model,
    ...     genome=genome
    ...     # cell_types=None means auto-generate from model output
    ... )
    """
    import crested
    
    # Handle input: file path or DataFrame
    if isinstance(bed_file, pd.DataFrame):
        bed_df = bed_file.copy()
        # Ensure standard column names if not already present
        if bed_df.shape[1] >= 3 and list(bed_df.columns[:3]) != ['chrom', 'start', 'end']:
            # Assume first 3 columns are chrom, start, end
            bed_df.columns = ["chrom", "start", "end"] + [f"col{i}" for i in range(3, bed_df.shape[1])]
    elif isinstance(bed_file, str):
        # Read BED file
        bed_df = pd.read_csv(bed_file, sep="\t", header=None)
    else:
        raise TypeError("bed_file must be either a file path (str) or a pandas DataFrame")
    
    # Parse BED format
    if bed_df.shape[1] < 3:
        raise ValueError("BED data must have at least 3 columns (chrom, start, end)")
    
    # Standardize column names if needed
    if 'chrom' not in bed_df.columns or 'start' not in bed_df.columns:
        bed_df.columns = ["chrom", "start", "end"] + [f"col{i}" for i in range(3, bed_df.shape[1])]
    
    # Check if ref/alt are provided
    has_ref_alt = bed_df.shape[1] >= 6
    if has_ref_alt:
        # Try to identify ref and alt columns
        # Common formats: chrom, start, end, name, ref, alt
        # Or: chrom, start, end, ref, alt
        if bed_df.shape[1] >= 6:
            # Assume columns 4 and 5 are ref/alt (0-indexed: 3 and 4 after chrom/start/end)
            bed_df["ref"] = bed_df.iloc[:, 4].astype(str).str.upper()
            bed_df["alt"] = bed_df.iloc[:, 5].astype(str).str.upper()
        elif bed_df.shape[1] >= 5:
            # Assume columns 3 and 4 are ref/alt
            bed_df["ref"] = bed_df.iloc[:, 3].astype(str).str.upper()
            bed_df["alt"] = bed_df.iloc[:, 4].astype(str).str.upper()
    
    # Infer seq_length from model if not provided
    if seq_length is None:
        try:
            seq_length = model.input_shape[1]
        except:
            raise ValueError("Could not infer seq_length from model. Please provide seq_length parameter.")
    
    import crested
    
    # Prepare all SNP data first (avoid per-SNP overhead)
    print(f"Preparing {len(bed_df)} SNP positions...")
    snp_data = []
    for idx, row in bed_df.iterrows():
        chrom = str(row["chrom"])
        snp_pos = int(row["start"])  # BED is 0-based
        
        # Calculate region centered on SNP
        center = snp_pos
        start = center - seq_length // 2
        end = start + seq_length
        
        # Determine which mutations to test
        has_ref = pd.notna(row.get("ref"))
        has_alt = pd.notna(row.get("alt"))
        
        snp_data.append({
            'chrom': chrom,
            'snp_pos': snp_pos,
            'start': start,
            'end': end,
            'ref': str(row["ref"]).upper() if has_ref else None,
            'alt': str(row["alt"]).upper() if has_alt else None,
        })
    
    # Process in chunks to avoid memory issues
    n_snps = len(snp_data)
    n_chunks = (n_snps + chunk_size - 1) // chunk_size
    print(f"Processing {n_snps} SNPs in {n_chunks} chunks of ~{chunk_size} SNPs each...")
    
    all_results = []
    
    for chunk_idx in range(n_chunks):
        start_idx = chunk_idx * chunk_size
        end_idx = min((chunk_idx + 1) * chunk_size, n_snps)
        chunk_snps = snp_data[start_idx:end_idx]
        
        print(f"Chunk {chunk_idx + 1}/{n_chunks}: Processing SNPs {start_idx}-{end_idx}...")
        
        # Extract sequences for this chunk
        sequences = []
        valid_snps = []
        for snp in chunk_snps:
            try:
                seq = genome.fetch(snp['chrom'], snp['start'], snp['end'])
                if len(seq) != seq_length:
                    continue
                sequences.append(seq)
                valid_snps.append(snp)
            except Exception:
                continue
        
        if not sequences:
            print(f"  No valid sequences in chunk {chunk_idx + 1}, skipping...")
            continue
        
        print(f"  Extracted {len(sequences)} valid sequences")
        
        # Batch predict all wildtype sequences for this chunk
        wt_preds = crested.tl.predict(input=sequences, model=model, batch_size=batch_size)
        wt_preds = np.asarray(wt_preds, dtype=np.float32)
        
        # Infer cell types from output if not provided (only need to do once)
        if chunk_idx == 0:
            if cell_types is None:
                n_cell_types = wt_preds.shape[1]
                cell_types_list = [f"CellType_{i}" for i in range(n_cell_types)]
            elif hasattr(cell_types, 'obs_names'):
                cell_types_list = list(map(str, cell_types.obs_names))
            else:
                cell_types_list = list(map(str, cell_types))
        
        # Build mutant sequences for this chunk
        all_mutants = []
        mutant_metadata = []
        
        for seq_idx, (seq, snp) in enumerate(zip(sequences, valid_snps)):
            snp_offset = snp['snp_pos'] - snp['start']
            if snp_offset < 0 or snp_offset >= len(seq):
                continue
                
            ref_from_seq = seq[snp_offset].upper()
            
            # Determine mutations to test
            if snp['ref'] and snp['alt']:
                ref_allele = snp['ref']
                if ref_allele != ref_from_seq:
                    ref_allele = ref_from_seq
                mutations = [(ref_allele, snp['alt'])]
            else:
                ref_allele = ref_from_seq
                bases = ["A", "C", "G", "T"]
                mutations = [(ref_allele, alt) for alt in bases if alt != ref_allele]
            
            if skip_ambiguous and ref_allele not in ["A", "C", "G", "T"]:
                continue
            
            for ref_base, alt_base in mutations:
                if alt_base == ref_base:
                    continue
                mutant_seq = seq[:snp_offset] + alt_base + seq[snp_offset+1:]
                all_mutants.append(mutant_seq)
                mutant_metadata.append({
                    'seq_idx': seq_idx,
                    'chrom': snp['chrom'],
                    'pos': snp['snp_pos'],
                    'ref': ref_base,
                    'alt': alt_base,
                })
        
        if not all_mutants:
            print(f"  No valid mutations in chunk {chunk_idx + 1}, skipping...")
            continue
        
        # Batch predict all mutants for this chunk
        print(f"  Predicting {len(all_mutants)} mutations...")
        mut_preds = crested.tl.predict(input=all_mutants, model=model, batch_size=batch_size)
        mut_preds = np.asarray(mut_preds, dtype=np.float32)
        
        # Build results for this chunk
        for mut_idx, (mut_pred, meta) in enumerate(zip(mut_preds, mutant_metadata)):
            seq_idx = meta['seq_idx']
            wt_pred = wt_preds[seq_idx]
            
            for ct_idx, ct_name in enumerate(cell_types_list):
                wt_val = wt_pred[ct_idx]
                mut_val = mut_pred[ct_idx]
                
                # Calculate log2fc
                epsilon = 1e-10
                log2fc = np.log2((mut_val + epsilon) / (wt_val + epsilon))
                delta = mut_val - wt_val
                
                mut_id = f"{meta['chrom']}:{meta['pos']}:{meta['ref']}>{meta['alt']}"
                
                all_results.append({
                    'chrom': meta['chrom'],
                    'pos': meta['pos'],
                    'ref': meta['ref'],
                    'alt': meta['alt'],
                    'mut_id': mut_id,
                    'cell_type': ct_name,
                    'wt_pred': wt_val,
                    'mut_pred': mut_val,
                    'log2fc': log2fc,
                    'delta': delta
                })
        
        print(f"  Chunk {chunk_idx + 1}/{n_chunks} complete: {len(all_results)} total results so far")
    
    # Combine all results
    if not all_results:
        return pd.DataFrame(columns=[
            "chrom", "pos", "ref", "alt", "mut_id", 
            "cell_type", "wt_pred", "mut_pred", "log2fc", "delta"
        ])
    
    print(f"Finished! Total results: {len(all_results)}")
    combined_df = pd.DataFrame(all_results)
    return combined_df