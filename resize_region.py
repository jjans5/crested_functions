def resize_region(region_str, length, summit=None):
    """
    Resize a region string like 'chr12:112924652-112925561' to `new_length`,
    centered on the middle of the region by default, or on `summit` if provided.
    """
    chrom, coords = region_str.split(':')
    start, end = map(int, coords.split('-'))

    center = summit if summit is not None else (start + end) // 2

    half_len = length // 2
    new_start = center - half_len
    new_end = center + half_len

    return f"{chrom}:{new_start}-{new_end}"
    
# similar example but with region names as input
regions_of_interest = [
    "chr21:18402602-18404716"
]  # FIRE enhancer region (Microglia enhancer)

regions_of_interest = [resize_region(x,length=2114) for x in regions_of_interest]
classes_of_interest = list(adata.obs_names)
class_idx = list(adata.obs_names.get_indexer(classes_of_interest))

scores, one_hot_encoded_sequences = crested.tl.contribution_scores(
    regions_of_interest,
    target_idx=class_idx,
    model=model, batch_size=24
)