import numpy as np
import os
import tempfile
import shutil
import logging
from datetime import datetime
from typing import Optional, Union
from tqdm import tqdm
import warnings
import anndata as ad
import crested

try:
    import psutil
    PSUTIL_AVAILABLE = True
except ImportError:
    PSUTIL_AVAILABLE = False
    warnings.warn("psutil not available. Memory monitoring will be disabled.")

def predict_chunked(
    adata, 
    model, 
    batch_size: int = 16,
    chunk_size: int = 1000,
    layer_name: str = "pred_model",
    memory_monitoring: bool = True,
    disk_based_saving: bool = False,
    temp_dir: Optional[str] = None,
    auto_adjust_chunks: bool = True,
    max_memory_percent: float = 80.0,
    verbose: bool = True
):
    """
    Comprehensive chunked prediction function for large AnnData objects.
    
    Parameters:
    -----------
    adata : ad.AnnData
        The AnnData object to make predictions on
    model : 
        The trained CREsted model
    batch_size : int, default=16
        Batch size for prediction (keep small for memory efficiency)
    chunk_size : int, default=1000
        Number of genomic regions to process per chunk
    layer_name : str, default="pred_model"
        Name for the predictions layer in adata
    memory_monitoring : bool, default=True
        Enable memory usage monitoring and warnings
    disk_based_saving : bool, default=False
        Save intermediate chunks to disk (for very large datasets)
    temp_dir : str, optional
        Directory for temporary files. If None and disk_based_saving=True, 
        creates a temporary directory with timestamp
    auto_adjust_chunks : bool, default=True
        Automatically adjust chunk_size based on available memory
    max_memory_percent : float, default=80.0
        Maximum memory usage percentage before reducing chunk size
    verbose : bool, default=True
        Enable detailed logging
        
    Returns:
    --------
    ad.AnnData
        The original adata with predictions added to layers
    """
    
    # Set up logging
    if verbose:
        log = logging.getLogger(__name__)
        if not log.handlers:
            logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    else:
        log = logging.getLogger(__name__)
        log.setLevel(logging.WARNING)
    
    # Always chunk by regions (vars)
    n_regions = adata.n_vars
    log.info(f"🚀 Starting chunked prediction for {adata.n_obs:,} cell types × {n_regions:,} regions")
    log.info(f"📊 Chunking by REGIONS - each chunk will have ALL {adata.n_obs:,} cell types")
    
    original_chunk_size = chunk_size
    log.info(f"📊 Initial settings: chunk_size={chunk_size:,}, batch_size={batch_size}")
    
    # Memory monitoring setup
    memory_before = None
    if memory_monitoring and PSUTIL_AVAILABLE:
        memory_before = psutil.virtual_memory()
        log.info(f"💾 Memory before prediction: {memory_before.percent:.1f}% used")
        log.info(f"💾 Available memory: {memory_before.available / (1024**3):.2f} GB")
        
        # Auto-adjust chunk size based on memory
        if auto_adjust_chunks:
            if memory_before.percent > max_memory_percent:
                chunk_size = max(chunk_size // 4, 100)
                log.warning(f"⚠️  High memory usage detected. Reducing chunk_size to {chunk_size:,}")
            elif memory_before.available < 4 * (1024**3):  # Less than 4GB available
                chunk_size = max(chunk_size // 2, 200)
                log.warning(f"⚠️  Limited memory available. Reducing chunk_size to {chunk_size:,}")
    
    elif memory_monitoring and not PSUTIL_AVAILABLE:
        log.warning("⚠️  Memory monitoring requested but psutil not available")
    
    # Setup temporary directory for disk-based saving
    temp_dir_created = False
    temp_dir_path = None
    
    if disk_based_saving:
        if temp_dir is None:
            # Create temp dir with timestamp to avoid conflicts
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")[:-3]  # microseconds to milliseconds
            temp_dir_name = f"crested_predictions_{timestamp}"
            temp_dir_path = tempfile.mkdtemp(prefix=temp_dir_name + "_")
            temp_dir_created = True
            log.info(f"📁 Created temporary directory: {temp_dir_path}")
        else:
            temp_dir_path = temp_dir
            os.makedirs(temp_dir_path, exist_ok=True)
        log.info(f"💾 Disk-based saving enabled. Using directory: {temp_dir_path}")
    
    # Calculate number of chunks
    n_chunks = (n_regions + chunk_size - 1) // chunk_size
    log.info(f"📦 Will process {n_chunks:,} chunks of size ≤{chunk_size:,} regions")
    
    # Storage for predictions
    all_predictions = []
    chunk_files = []
    
    # Progress bar setup
    chunk_pbar = tqdm(
        total=n_chunks, 
        desc="Processing regions", 
        unit="chunk",
        bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} chunks [{elapsed}<{remaining}, {rate_fmt}]"
    )
    
    try:
        for i in range(n_chunks):
            start_idx = i * chunk_size
            end_idx = min((i + 1) * chunk_size, n_regions)
            chunk_size_actual = end_idx - start_idx
            
            chunk_pbar.set_description(f"Chunk {i+1}/{n_chunks} ({chunk_size_actual:,} regions)")
            
            if verbose:
                log.info(f"📊 Processing chunk {i+1}/{n_chunks}: regions {start_idx:,}-{end_idx:,} ({chunk_size_actual:,} regions)")
            
            # Memory check before processing chunk
            if memory_monitoring and PSUTIL_AVAILABLE:
                current_memory = psutil.virtual_memory()
                if current_memory.percent > max_memory_percent:
                    log.warning(f"⚠️  Memory usage high ({current_memory.percent:.1f}%). Consider reducing chunk_size further.")
                if verbose:
                    chunk_pbar.set_postfix({"Memory": f"{current_memory.percent:.1f}%"})
            
            # Extract chunk - ALL cell types, subset of regions
            adata_chunk = adata[:, start_idx:end_idx].copy()
            
            try:
                # Make predictions on chunk
                if verbose:
                    log.info(f"  🔮 Making predictions on {adata_chunk.n_obs:,} obs × {adata_chunk.n_vars:,} vars...")
                
                predictions_chunk = crested.tl.predict(
                    adata_chunk, 
                    model, 
                    batch_size=batch_size
                )
                
                # Store predictions in the chunk's layers
                adata_chunk.layers[layer_name] = predictions_chunk.T  # Store as (n_vars, n_obs)
                
                if verbose:
                    log.info(f"  ✅ Predictions added to chunk layers with shape: {adata_chunk.layers[layer_name].shape}")
                
                # Handle storage based on options
                if disk_based_saving:
                    # Save chunk as h5ad file
                    chunk_file = os.path.join(temp_dir_path, f"chunk_{i:05d}_predictions.h5ad")
                    adata_chunk.write(chunk_file)
                    chunk_files.append(chunk_file)
                    if verbose:
                        log.info(f"  💾 Saved chunk to: {os.path.basename(chunk_file)}")
                else:
                    # Keep in memory
                    all_predictions.append(adata_chunk)
                
            except Exception as e:
                log.error(f"❌ Error processing chunk {i+1}: {e}")
                raise
            finally:
                # Clean up chunk only if using disk-based saving
                if disk_based_saving:
                    del adata_chunk
            
            # Update progress bar
            chunk_pbar.update(1)
            
            # Memory status update for progress bar
            if memory_monitoring and PSUTIL_AVAILABLE:
                current_memory = psutil.virtual_memory()
                chunk_pbar.set_postfix({"Memory": f"{current_memory.percent:.1f}%"})
        
        # Close progress bar
        chunk_pbar.close()
        
        # Concatenate all predictions using AnnData
        log.info("🔗 Concatenating all predictions using AnnData...")
        
        if disk_based_saving:
            # Load AnnData chunks from disk
            log.info(f"📂 Loading {len(chunk_files):,} AnnData files from disk...")
            
            # Progress bar for loading files
            load_pbar = tqdm(
                chunk_files, 
                desc="Loading chunks", 
                unit="file",
                bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} files [{elapsed}<{remaining}]"
            )
            
            all_predictions = []
            for chunk_file in load_pbar:
                chunk_adata = ad.read_h5ad(chunk_file)
                all_predictions.append(chunk_adata)
                load_pbar.set_postfix({"Shape": f"{chunk_adata.shape}"})
            
            load_pbar.close()
        
        # Concatenate all chunks along vars (axis=1) using AnnData
        log.info("🔗 Performing AnnData concatenation along vars...")
        adata_concat = ad.concat(all_predictions, axis=1, merge="same")
        
        log.info(f"✅ Concatenated AnnData shape: {adata_concat.shape}")
        log.info(f"🎯 Expected shape: {adata.shape}")
        
        # Validate shapes
        if adata_concat.shape != adata.shape:
            raise ValueError(f"❌ Concatenated shape mismatch: got {adata_concat.shape} expected {adata.shape}")
        
        # Copy predictions to original adata
        log.info(f"📊 Adding predictions to adata.layers['{layer_name}']...")
        adata.layers[layer_name] = adata_concat.layers[layer_name]
        
        log.info(f"🎉 Successfully added predictions to adata.layers['{layer_name}']")
        
        # Final memory report
        if memory_monitoring and PSUTIL_AVAILABLE and memory_before:
            memory_after = psutil.virtual_memory()
            log.info(f"💾 Memory after prediction: {memory_after.percent:.1f}% used")
            memory_change = memory_after.percent - memory_before.percent
            log.info(f"💾 Memory change: {memory_change:+.1f}%")
        
        return adata
        
    except Exception as e:
        log.error(f"❌ Prediction failed: {e}")
        raise
        
    finally:
        # Ensure progress bar is closed
        if 'chunk_pbar' in locals():
            chunk_pbar.close()
        if 'load_pbar' in locals():
            load_pbar.close()
        
        # Cleanup temporary files and directory
        if disk_based_saving and chunk_files:
            log.info(f"🧹 Cleaning up {len(chunk_files):,} temporary files...")
            cleanup_pbar = tqdm(
                chunk_files, 
                desc="Cleaning files", 
                unit="file",
                leave=False
            )
            
            for chunk_file in cleanup_pbar:
                try:
                    if os.path.exists(chunk_file):
                        os.remove(chunk_file)
                except OSError as e:
                    log.warning(f"⚠️  Could not remove temporary file {chunk_file}: {e}")
            
            cleanup_pbar.close()
        
        # Cleanup temporary directory if we created it
        if disk_based_saving and temp_dir_created and temp_dir_path:
            try:
                if os.path.exists(temp_dir_path):
                    shutil.rmtree(temp_dir_path)
                    log.info(f"🧹 Cleaned up temporary directory: {os.path.basename(temp_dir_path)}")
            except OSError as e:
                log.warning(f"⚠️  Could not clean up temporary directory {temp_dir_path}: {e}")
