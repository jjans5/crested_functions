# Your input dictionary
regions_of_interest = {
    'OAS1': 'chr12:112924652-112925561',
    'MX1': 'chr21:41411243-41412108',
    'IL32': 'chr16:3103284-3104174',
    'IL18': 'chr11:112157030-112157939',
    'WSCD1_pk4': 'chr17:6110198-6111245',
    'WSCD1_pk3': 'chr17:6108525-6108695',
    'WSCD1_pk2': 'chr17:6104633-6105474',
    'WSCD1_pk1': 'chr17:6101838-6103814',
    'WSCD1_pk0': 'chr17:6090124-6090753'
    
}

# Generate new resized region dictionary
resized_regions = {
    gene: resize_region(region, length=2114)
    for gene, region in regions_of_interest.items()
}

# If you want to inspect it
for gene, resized in resized_regions.items():
    print(f"{gene}: {resized}")

import os
import re
import pandas as pd

outdir = "ism_results"
os.makedirs(outdir, exist_ok=True)

region_re = re.compile(r'^(?P<chrom>[^:]+):(?P<start>\d+)-(?P<end>\d+)$')

for enhancer_name, region in resized_regions.items():
    m = region_re.match(region)
    if not m:
        print(f"[skip] Could not parse region '{region}' for {enhancer_name}")
        continue

    chrom = m.group("chrom")
    start = int(m.group("start"))
    end   = int(m.group("end"))  # pysam FastaFile.fetch uses end-exclusive

    # Fetch sequence and run ISM
    seq = genome.fetch(chrom, start, end)
    res = insilico_mutagenesis_vect(
        seq=seq,
        model=model,
        adata=adata,
        chrom=chrom,
        start=start,           # genomic coordinate of seq[0]
        skip_ambiguous=True,
        return_long=True
    )

    wt_pred     = res["wt_pred"]     # pd.Series
    mut_df_wide = res["mut_df_wide"] # DataFrame (wide)
    mut_df_long = res["mut_df_long"] # DataFrame (long + delta)

    # Strongest negative-impact mutations per cell type (delta ascending)
    top_loss = (mut_df_long.sort_values(["cell_type", "delta"])
                          .groupby("cell_type", group_keys=False)
                          .head(20))

    # --- Save as TSVs ---
    # 1) WT predictions: two columns [cell_type, value]
    wt_path = os.path.join(outdir, f"{enhancer_name}__wt.tsv")
    wt_pred.rename("value").to_frame().to_csv(wt_path, sep="\t")

    # 2) Mutant predictions (wide)
    wide_path = os.path.join(outdir, f"{enhancer_name}__mut_wide.tsv")
    mut_df_wide.to_csv(wide_path, sep="\t", index=False)

    # 3) Mutant predictions (long + delta)
    long_path = os.path.join(outdir, f"{enhancer_name}__mut_long.tsv")
    mut_df_long.to_csv(long_path, sep="\t", index=False)

    # 4) Top losses per cell type
    toploss_path = os.path.join(outdir, f"{enhancer_name}__top_loss.tsv")
    top_loss.to_csv(toploss_path, sep="\t", index=False)

    print(f"[done] {enhancer_name} → {wt_path}, {wide_path}, {long_path}, {toploss_path}")