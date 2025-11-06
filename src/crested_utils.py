"""
Utility functions for working with CREsted predictions across species.

This module provides functions to:
1. Make predictions on AnnData objects (with chunking support for large datasets)
2. Align cell types between species for comparison
3. Calculate correlations between predicted and observed values
4. Resize genomic regions
"""

import numpy as np
import pandas as pd
import anndata as ad
import crested
import warnings
import logging
from typing import Optional, Union, List, Tuple
from tqdm import tqdm
import os
import tempfile
import shutil
from datetime import datetime

try:
    import psutil
    PSUTIL_AVAILABLE = True
except ImportError:
    PSUTIL_AVAILABLE = False
    warnings.warn("psutil not available. Memory monitoring will be disabled.")


def predict_regions(
    model,
    regions: Optional[Union[List[str], ad.AnnData]] = None,
    adata: Optional[ad.AnnData] = None,
    genome: Optional[crested.Genome] = None,
    batch_size: int = 16,
    chunk_size: int = 1000,
    layer_name: str = "predicted",
    target_cell_types: Optional[List[str]] = None,
    align_cell_types: bool = True,
    fill_missing_with_zeros: bool = True,
    disk_based_saving: bool = False,
    temp_dir: Optional[str] = None,
    verbose: bool = True,
) -> ad.AnnData:
    """
    Make predictions on genomic regions and return an AnnData object.
    
    This function handles predictions for:
    - A list of region strings (e.g., ['chr1:1000-2000', ...])
    - An existing AnnData object with regions in .var_names
    
    Parameters
    ----------
    model : CREsted model
        The trained CREsted model for making predictions
    regions : list of str or AnnData, optional
        Either a list of region strings or an AnnData object.
        If AnnData, will use its regions and cell type information.
    adata : AnnData, optional
        Alternative way to pass AnnData object (deprecated, use regions instead)
    genome : crested.Genome, optional
        Genome object required if regions is a list of strings
    batch_size : int, default=16
        Batch size for prediction
    chunk_size : int, default=1000
        Number of regions to process per chunk (for memory efficiency)
    layer_name : str, default='predicted'
        Name for the predictions layer in output AnnData
    target_cell_types : list of str, optional
        Desired cell types in output (in specific order).
        If provided and align_cell_types=True, will align predictions to this list.
    align_cell_types : bool, default=True
        Whether to align cell types to target_cell_types order
    fill_missing_with_zeros : bool, default=True
        If True, missing cell types will be filled with zeros.
        If False, only common cell types will be retained.
    disk_based_saving : bool, default=False
        Save intermediate chunks to disk for very large datasets
    temp_dir : str, optional
        Directory for temporary files
    verbose : bool, default=True
        Enable detailed logging
        
    Returns
    -------
    AnnData
        AnnData object with:
        - .X: predicted values (n_cell_types x n_regions)
        - .layers[layer_name]: same as .X
        - .obs: cell type information
        - .var: region information
        
    Examples
    --------
    >>> # From region strings
    >>> regions = ['chr1:1000-2000', 'chr1:3000-4000']
    >>> adata_pred = predict_regions(model, regions=regions, genome=genome)
    
    >>> # From existing AnnData
    >>> adata_pred = predict_regions(model, regions=adata_species)
    
    >>> # With cell type alignment
    >>> target_cts = ['T_cell', 'B_cell', 'Macrophage']
    >>> adata_pred = predict_regions(
    ...     model, regions=adata_species,
    ...     target_cell_types=target_cts,
    ...     align_cell_types=True
    ... )
    """
    # Set up logging
    log = logging.getLogger(__name__)
    if verbose:
        if not log.handlers:
            logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    else:
        log.setLevel(logging.WARNING)
    
    # Handle input
    if adata is not None and regions is None:
        regions = adata
        warnings.warn(
            "Passing adata as a keyword argument is deprecated. "
            "Use regions=adata instead.",
            DeprecationWarning
        )
    
    if regions is None:
        raise ValueError("Must provide either 'regions' (list or AnnData)")
    
    # Case 1: regions is already an AnnData object
    if isinstance(regions, ad.AnnData):
        log.info(f"📊 Input is AnnData with shape {regions.shape}")
        adata_input = regions
        
    # Case 2: regions is a list of region strings
    elif isinstance(regions, (list, tuple)):
        log.info(f"📊 Processing {len(regions)} region strings")
        
        # Make predictions directly on the region strings
        # crested.tl.predict can take a list of region strings
        log.info("� Making predictions on regions...")
        predictions = crested.tl.predict(
            input=regions,
            model=model,
            batch_size=batch_size,
            genome=genome
        )
        
        # predictions shape: (n_regions, n_cell_types)
        predictions = np.asarray(predictions, dtype=np.float32)
        n_regions_pred, n_cell_types_pred = predictions.shape
        
        log.info(f"✅ Predictions shape: {predictions.shape} (regions x cell_types)")
        
        # Determine cell types
        if target_cell_types is not None:
            if len(target_cell_types) != n_cell_types_pred:
                log.warning(
                    f"⚠️ target_cell_types length ({len(target_cell_types)}) doesn't match "
                    f"prediction output ({n_cell_types_pred}). Using predicted cell types."
                )
                cell_types = [f"CellType_{i}" for i in range(n_cell_types_pred)]
            else:
                cell_types = target_cell_types
        else:
            # Auto-generate cell type names
            cell_types = [f"CellType_{i}" for i in range(n_cell_types_pred)]
            log.info(f"📝 Auto-generated {len(cell_types)} cell type names")
        
        # Create AnnData with predictions
        # Note: AnnData expects (n_obs x n_vars) = (n_cell_types x n_regions)
        var = pd.DataFrame(index=regions)
        obs = pd.DataFrame(index=cell_types)
        
        # Transpose predictions to (n_cell_types, n_regions)
        X = predictions.T
        adata_input = ad.AnnData(X=X, obs=obs, var=var)
        adata_input.layers[layer_name] = X
        
        # Return directly without going through _predict_chunked
        # since we already have the predictions
        if align_cell_types and target_cell_types is not None:
            log.info(f"🔄 Aligning cell types to target list of {len(target_cell_types)} cell types")
            adata_input = align_adata_cell_types(
                adata_input,
                target_cell_types=target_cell_types,
                fill_missing=fill_missing_with_zeros,
                verbose=verbose
            )
        
        return adata_input
        
    else:
        raise TypeError("regions must be either a list of strings or an AnnData object")
    
    # Make predictions using chunked approach
    adata_with_pred = _predict_chunked(
        adata=adata_input,
        model=model,
        batch_size=batch_size,
        chunk_size=chunk_size,
        layer_name=layer_name,
        disk_based_saving=disk_based_saving,
        temp_dir=temp_dir,
        verbose=verbose
    )
    
    # Handle cell type alignment if requested
    if align_cell_types and target_cell_types is not None:
        log.info(f"🔄 Aligning cell types to target list of {len(target_cell_types)} cell types")
        adata_with_pred = align_adata_cell_types(
            adata_with_pred,
            target_cell_types=target_cell_types,
            fill_missing=fill_missing_with_zeros,
            verbose=verbose
        )
    
    return adata_with_pred


def align_adata_cell_types(
    adata: ad.AnnData,
    target_cell_types: List[str],
    fill_missing: bool = True,
    verbose: bool = True
) -> ad.AnnData:
    """
    Align AnnData cell types (obs) to a target list of cell types.
    
    This is useful when comparing predictions across species where not all
    cell types may be present in each species.
    
    Parameters
    ----------
    adata : AnnData
        Input AnnData object
    target_cell_types : list of str
        Desired cell types in specific order
    fill_missing : bool, default=True
        If True, missing cell types will be added with zeros.
        If False, only common cell types will be retained.
    verbose : bool, default=True
        Enable logging
        
    Returns
    -------
    AnnData
        New AnnData with aligned cell types
        
    Examples
    --------
    >>> target_cts = ['T_cell', 'B_cell', 'Macrophage']
    >>> adata_aligned = align_adata_cell_types(adata_species, target_cts)
    """
    log = logging.getLogger(__name__)
    if verbose:
        log.setLevel(logging.INFO)
    
    # Identify present and missing cell types
    present = [ct for ct in target_cell_types if ct in adata.obs_names]
    missing = [ct for ct in target_cell_types if ct not in adata.obs_names]
    
    log.info(f"📊 Cell types: {len(present)} present, {len(missing)} missing")
    
    if not fill_missing and missing:
        # Only keep common cell types
        log.info(f"⚠️  Removing {len(missing)} missing cell types")
        return adata[present, :].copy()
    
    # Subset to present cell types
    adata_present = adata[present, :].copy()
    
    # Keep X as dense float32
    X_present = np.asarray(adata_present.X)
    if X_present.dtype != np.float32:
        X_present = X_present.astype(np.float32, copy=False)
    
    # Build zeros for missing cell types
    if missing:
        log.info(f"➕ Adding {len(missing)} missing cell types with zeros")
        X_missing = np.zeros((len(missing), adata.n_vars), dtype=X_present.dtype)
        X_new = np.vstack([X_present, X_missing])
    else:
        X_new = X_present
    
    # Build obs DataFrame in present+missing order
    obs_temp = pd.concat([
        adata_present.obs,
        pd.DataFrame(index=pd.Index(missing, name=adata.obs_names.name)),
    ], axis=0)
    
    # Copy var
    var_new = adata.var.copy()
    
    # Create temporary AnnData
    adata_temp = ad.AnnData(X=X_new, obs=obs_temp, var=var_new)
    
    # Preserve layers
    for layer_key in adata_present.layers.keys():
        L = np.asarray(adata_present.layers[layer_key])
        if missing:
            L_missing = np.zeros((len(missing), adata.n_vars), dtype=L.dtype)
            L_new = np.vstack([L, L_missing])
        else:
            L_new = L
        adata_temp.layers[layer_key] = L_new
    
    # Preserve obsm (e.g., weights)
    for obsm_key in adata_present.obsm.keys():
        W = np.asarray(adata_present.obsm[obsm_key])
        if W.ndim == 1:
            W = W.reshape(-1, 1)
        if missing:
            W_new = np.vstack([W, np.zeros((len(missing), W.shape[1]), dtype=W.dtype)])
        else:
            W_new = W
        adata_temp.obsm[obsm_key] = W_new
    
    # Reorder to exact target cell type order
    adata_aligned = adata_temp[target_cell_types, :].copy()
    
    assert adata_aligned.X.dtype == np.float32
    log.info(f"✅ Aligned AnnData shape: {adata_aligned.shape}")
    
    return adata_aligned


def _predict_chunked(
    adata: ad.AnnData,
    model,
    batch_size: int = 16,
    chunk_size: int = 1000,
    layer_name: str = "predicted",
    memory_monitoring: bool = True,
    disk_based_saving: bool = False,
    temp_dir: Optional[str] = None,
    auto_adjust_chunks: bool = True,
    max_memory_percent: float = 80.0,
    verbose: bool = True
) -> ad.AnnData:
    """
    Internal function for chunked prediction on AnnData.
    
    Processes large AnnData objects in chunks to avoid memory issues.
    Chunks are created along the regions (vars) dimension.
    """
    log = logging.getLogger(__name__)
    
    n_regions = adata.n_vars
    log.info(f"🚀 Starting chunked prediction for {adata.n_obs:,} cell types × {n_regions:,} regions")
    log.info(f"📊 Initial settings: chunk_size={chunk_size:,}, batch_size={batch_size}")
    
    # Memory monitoring
    memory_before = None
    if memory_monitoring and PSUTIL_AVAILABLE:
        memory_before = psutil.virtual_memory()
        log.info(f"💾 Memory before prediction: {memory_before.percent:.1f}% used")
        log.info(f"💾 Available memory: {memory_before.available / (1024**3):.2f} GB")
        
        if auto_adjust_chunks:
            if memory_before.percent > max_memory_percent:
                chunk_size = max(chunk_size // 4, 100)
                log.warning(f"⚠️  High memory usage. Reducing chunk_size to {chunk_size:,}")
            elif memory_before.available < 4 * (1024**3):
                chunk_size = max(chunk_size // 2, 200)
                log.warning(f"⚠️  Limited memory. Reducing chunk_size to {chunk_size:,}")
    
    # Setup temporary directory
    temp_dir_created = False
    temp_dir_path = None
    
    if disk_based_saving:
        if temp_dir is None:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")[:-3]
            temp_dir_name = f"crested_predictions_{timestamp}"
            temp_dir_path = tempfile.mkdtemp(prefix=temp_dir_name + "_")
            temp_dir_created = True
            log.info(f"📁 Created temporary directory: {temp_dir_path}")
        else:
            temp_dir_path = temp_dir
            os.makedirs(temp_dir_path, exist_ok=True)
    
    # Calculate chunks
    n_chunks = (n_regions + chunk_size - 1) // chunk_size
    log.info(f"📦 Will process {n_chunks:,} chunks of size ≤{chunk_size:,} regions")
    
    # Storage
    all_predictions = []
    chunk_files = []
    
    # Progress bar
    chunk_pbar = tqdm(
        total=n_chunks,
        desc="Processing regions",
        unit="chunk",
        bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} chunks [{elapsed}<{remaining}]"
    )
    
    try:
        for i in range(n_chunks):
            start_idx = i * chunk_size
            end_idx = min((i + 1) * chunk_size, n_regions)
            chunk_size_actual = end_idx - start_idx
            
            chunk_pbar.set_description(f"Chunk {i+1}/{n_chunks} ({chunk_size_actual:,} regions)")
            
            # Memory check
            if memory_monitoring and PSUTIL_AVAILABLE:
                current_memory = psutil.virtual_memory()
                if current_memory.percent > max_memory_percent:
                    log.warning(f"⚠️  Memory usage high ({current_memory.percent:.1f}%)")
                chunk_pbar.set_postfix({"Memory": f"{current_memory.percent:.1f}%"})
            
            # Extract chunk (all cell types, subset of regions)
            adata_chunk = adata[:, start_idx:end_idx].copy()
            
            try:
                # Make predictions
                predictions_chunk = crested.tl.predict(
                    adata_chunk,
                    model,
                    batch_size=batch_size
                )
                
                # Store predictions - note: crested returns (n_regions, n_celltypes)
                # We need (n_celltypes, n_regions) for consistency
                adata_chunk.layers[layer_name] = predictions_chunk.T
                
                # Handle storage
                if disk_based_saving:
                    chunk_file = os.path.join(temp_dir_path, f"chunk_{i:05d}.h5ad")
                    adata_chunk.write(chunk_file)
                    chunk_files.append(chunk_file)
                else:
                    all_predictions.append(adata_chunk)
                
            except Exception as e:
                log.error(f"❌ Error processing chunk {i+1}: {e}")
                raise
            finally:
                if disk_based_saving:
                    del adata_chunk
            
            chunk_pbar.update(1)
        
        chunk_pbar.close()
        
        # Concatenate predictions
        log.info("🔗 Concatenating all predictions...")
        
        if disk_based_saving:
            log.info(f"📂 Loading {len(chunk_files):,} files from disk...")
            all_predictions = [ad.read_h5ad(f) for f in tqdm(chunk_files, desc="Loading")]
        
        # Concatenate along vars axis
        adata_concat = ad.concat(all_predictions, axis=1, merge="same")
        log.info(f"✅ Concatenated shape: {adata_concat.shape}")
        
        # Validate
        if adata_concat.shape != adata.shape:
            raise ValueError(f"❌ Shape mismatch: got {adata_concat.shape} expected {adata.shape}")
        
        # Copy predictions to original adata
        adata.layers[layer_name] = adata_concat.layers[layer_name]
        
        # Final memory report
        if memory_monitoring and PSUTIL_AVAILABLE and memory_before:
            memory_after = psutil.virtual_memory()
            log.info(f"💾 Memory after: {memory_after.percent:.1f}% used")
            log.info(f"💾 Memory change: {memory_after.percent - memory_before.percent:+.1f}%")
        
        log.info(f"🎉 Successfully added predictions to adata.layers['{layer_name}']")
        return adata
        
    except Exception as e:
        log.error(f"❌ Prediction failed: {e}")
        raise
        
    finally:
        # Cleanup
        if 'chunk_pbar' in locals():
            chunk_pbar.close()
        
        if disk_based_saving and chunk_files:
            log.info(f"🧹 Cleaning up {len(chunk_files):,} temporary files...")
            for chunk_file in chunk_files:
                try:
                    if os.path.exists(chunk_file):
                        os.remove(chunk_file)
                except OSError as e:
                    log.warning(f"⚠️  Could not remove {chunk_file}: {e}")
        
        if disk_based_saving and temp_dir_created and temp_dir_path:
            try:
                if os.path.exists(temp_dir_path):
                    shutil.rmtree(temp_dir_path)
                    log.info(f"🧹 Cleaned up temporary directory")
            except OSError as e:
                log.warning(f"⚠️  Could not clean up {temp_dir_path}: {e}")


def rowwise_correlation(
    A: pd.DataFrame,
    B: pd.DataFrame,
    method: str = "pearson"
) -> pd.DataFrame:
    """
    Compute row-wise correlation between two DataFrames.
    
    Returns an (n x n) correlation matrix R where R[i, j] is the correlation
    between row i of A and row j of B, computed across their common columns.
    
    Parameters
    ----------
    A : DataFrame
        First DataFrame (n_rows_A x n_columns)
    B : DataFrame
        Second DataFrame (n_rows_B x n_columns)
    method : str, default='pearson'
        Correlation method: 'pearson' or 'spearman'
        
    Returns
    -------
    DataFrame
        Correlation matrix (n_rows_A x n_rows_B)
        
    Examples
    --------
    >>> # Compare predicted vs true for common cell types
    >>> pred = pd.DataFrame(adata.layers['predicted'], index=adata.obs_names)
    >>> true = pd.DataFrame(adata.X, index=adata.obs_names)
    >>> corr = rowwise_correlation(true, pred)
    >>> # Diagonal = self-correlation (how well each cell type is predicted)
    >>> print(corr.diagonal())
    """
    # Align on common columns
    cols = A.columns.intersection(B.columns)
    if len(cols) == 0:
        raise ValueError("A and B have no columns in common")
    
    A_ = A[cols]
    B_ = B[cols]
    
    # Spearman: rank transform
    if method.lower() == "spearman":
        A_ = A_.rank(axis=1, method="average")
        B_ = B_.rank(axis=1, method="average")
    elif method.lower() != "pearson":
        raise ValueError("method must be 'pearson' or 'spearman'")
    
    # Convert to numpy
    X = A_.to_numpy(dtype=float)
    Y = B_.to_numpy(dtype=float)
    m = X.shape[1]
    
    # Center rows
    X0 = X - X.mean(axis=1, keepdims=True)
    Y0 = Y - Y.mean(axis=1, keepdims=True)
    
    # Std per row (protect against zeros)
    Xstd = X0.std(axis=1, ddof=1, keepdims=True)
    Ystd = Y0.std(axis=1, ddof=1, keepdims=True)
    Xstd[Xstd == 0] = 1.0
    Ystd[Ystd == 0] = 1.0
    
    # Z-score
    Xz = X0 / Xstd
    Yz = Y0 / Ystd
    
    # Correlation via matrix multiply
    R = (Xz @ Yz.T) / (m - 1)
    
    return pd.DataFrame(R, index=A.index, columns=B.index)


def resize_region(
    region: str,
    new_length: int,
    summit: Optional[int] = None
) -> str:
    """
    Resize a genomic region to a new length.
    
    Parameters
    ----------
    region : str
        Region string in format 'chr:start-end' (e.g., 'chr1:1000-2000')
    new_length : int
        Desired length of the resized region
    summit : int, optional
        Position to center the region on. If None, uses the midpoint.
        
    Returns
    -------
    str
        Resized region string
        
    Examples
    --------
    >>> resize_region('chr1:1000-2000', 500)
    'chr1:1250-1750'
    >>> resize_region('chr1:1000-2000', 500, summit=1600)
    'chr1:1350-1850'
    """
    chrom, coords = region.split(':')
    start, end = map(int, coords.split('-'))
    
    center = summit if summit is not None else (start + end) // 2
    
    half_len = new_length // 2
    new_start = center - half_len
    new_end = center + half_len
    
    return f"{chrom}:{new_start}-{new_end}"


def compare_predictions(
    adata: ad.AnnData,
    prediction_layer: str = "predicted",
    true_layer: Optional[str] = None,
    cell_types: Optional[List[str]] = None,
    method: str = "pearson"
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Compare predicted vs true values for an AnnData object.
    
    Parameters
    ----------
    adata : AnnData
        AnnData with predictions in layers
    prediction_layer : str, default='predicted'
        Name of layer containing predictions
    true_layer : str, optional
        Name of layer containing true values. If None, uses .X
    cell_types : list of str, optional
        Cell types to compare. If None, uses all non-zero cell types.
    method : str, default='pearson'
        Correlation method: 'pearson' or 'spearman'
        
    Returns
    -------
    corr_df : DataFrame
        Row-wise correlation matrix
    diag_corr : DataFrame
        Diagonal correlations (self-correlation per cell type)
        
    Examples
    --------
    >>> corr_matrix, self_corr = compare_predictions(adata_with_pred)
    >>> print(f"Mean self-correlation: {self_corr['correlation'].mean():.3f}")
    """
    # Get true values
    if true_layer is None:
        true_data = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names)
    else:
        true_data = pd.DataFrame(adata.layers[true_layer], index=adata.obs_names, columns=adata.var_names)
    
    # Get predictions
    pred_data = pd.DataFrame(adata.layers[prediction_layer], index=adata.obs_names, columns=adata.var_names)
    
    # Filter to specified cell types or non-zero cell types
    if cell_types is None:
        # Keep cell types with non-zero values
        non_zero_mask = (true_data.sum(axis=1) > 0) | (pred_data.sum(axis=1) > 0)
        cell_types = adata.obs_names[non_zero_mask].tolist()
    
    # Subset
    true_data_sub = true_data.loc[cell_types]
    pred_data_sub = pred_data.loc[cell_types]
    
    # Compute correlations
    corr_df = rowwise_correlation(true_data_sub, pred_data_sub, method=method)
    
    # Extract diagonal (self-correlations)
    diag_values = np.diag(corr_df.values)
    diag_corr = pd.DataFrame({
        'cell_type': cell_types,
        'correlation': diag_values
    })
    
    return corr_df, diag_corr


# Backward compatibility aliases
predict_chunked = _predict_chunked
rowwise_corr = rowwise_correlation
