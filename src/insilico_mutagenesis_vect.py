import numpy as np
import pandas as pd

def insilico_mutagenesis_vect(
    seq: str,
    model,
    adata,
    chrom: str | None = None,
    start: int | None = None,   # genomic start for seq[0]; if None we report 0-based pos0
    skip_ambiguous: bool = True,
    batch_size=2,
    return_long: bool = True,   # also return tidy long table with deltas
):
    """
    Vectorized single-nucleotide in-silico mutagenesis (ISM).

    Assumptions
    -----------
    - crested.tl.predict(input=[str, ...], model=model) -> np.ndarray (N x C)
    - Columns correspond to `adata.obs_names` (length C)

    Returns
    -------
    dict with:
      wt_pred:      pd.Series (index=adata.obs_names)
      mut_df_wide:  pd.DataFrame with metadata columns ['chrom','genomic_pos','pos0','ref','alt','mut_id']
                    plus one column per cell type (absolute predictions)
      mut_df_long:  (if return_long) tidy DataFrame with columns
                    ['chrom','genomic_pos','pos0','ref','alt','mut_id','cell_type','value','delta']
    """
    import crested

    seq = seq.upper()
    bases = ("A", "C", "G", "T")
    cell_types = list(map(str, adata.obs_names))
    L = len(seq)
    C = len(cell_types)

    # 1) WT prediction
    wt = crested.tl.predict(input=[seq], model=model,batch_size=batch_size)
    wt = np.asarray(wt, dtype=np.float32).reshape(1, -1)
    if wt.shape[1] != C:
        raise ValueError(f"Prediction size mismatch: got {wt.shape[1]} cell types, expected {C} from adata.obs_names")
    wt_pred = pd.Series(wt[0], index=cell_types, name="wildtype")

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
        mut_df_wide = pd.DataFrame(columns=["chrom","genomic_pos","pos0","ref","alt","mut_id"] + cell_types)
        out = {"wt_pred": wt_pred, "mut_df_wide": mut_df_wide}
        if return_long:
            out["mut_df_long"] = pd.DataFrame(columns=["chrom","genomic_pos","pos0","ref","alt","mut_id","cell_type","value","delta"])
        return out

    # 3) Predict all mutants in a single call
    preds = crested.tl.predict(input=mutant_seqs, model=model,batch_size=batch_size)
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

    pred_df = pd.DataFrame(preds, columns=cell_types)
    mut_df_wide = pd.concat([meta_df, pred_df], axis=1)

    out = {"wt_pred": wt_pred, "mut_df_wide": mut_df_wide}

    # 5) Optional: tidy long table with deltas
    if return_long:
        long = mut_df_wide.melt(
            id_vars=["chrom","genomic_pos","pos0","ref","alt","mut_id"],
            value_vars=cell_types,
            var_name="cell_type",
            value_name="value"
        )
        # map wt once (fast)
        wt_map = wt_pred.to_dict()
        long["delta"] = (long["value"] - long["cell_type"].map(wt_map)).astype(np.float32)
        out["mut_df_long"] = long

    return out