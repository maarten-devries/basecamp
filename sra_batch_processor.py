#!/usr/bin/env python3
"""
Module for processing large batches of SRA IDs efficiently.
This module provides functions that can be imported and used in a Jupyter notebook
to convert SRA IDs to BioProject IDs and GEO/ArrayExpress identifiers.
"""

import os
import sys
import logging
import time
from typing import List, Dict, Optional, Union
import pandas as pd
from pathlib import Path
import importlib.util

# Try to determine if we're in a Jupyter notebook
def is_notebook():
    try:
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':  # Jupyter notebook or qtconsole
            return True
        elif shell == 'TerminalInteractiveShell':  # Terminal IPython
            return False
        else:
            return False  # Other type
    except NameError:
        return False  # Probably standard Python interpreter

# Import appropriate tqdm based on environment
if is_notebook():
    try:
        from tqdm.notebook import tqdm
    except ImportError:
        from tqdm import tqdm
else:
    from tqdm import tqdm

from sra_id_converter import convert_sra_ids

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def read_sra_ids_from_file(file_path: str) -> List[str]:
    """
    Read SRA IDs from a file, one ID per line.
    
    Args:
        file_path: Path to the file containing SRA IDs
        
    Returns:
        List of SRA IDs
    """
    with open(file_path, 'r') as f:
        # Strip whitespace and filter out empty lines
        return [line.strip() for line in f if line.strip()]

def process_sra_ids(
    sra_ids: List[str],
    batch_size: int = 10,
    max_workers: int = 2,
    delay_between_batches: float = 3.0,
    cache_file: Optional[str] = "sra_cache.json",
    max_retries: int = 3,
    retry_delay: float = 5.0,
    show_progress: bool = True
) -> Dict[str, Dict[str, str]]:
    """
    Process a large list of SRA IDs efficiently with batching and caching.
    
    Args:
        sra_ids: List of SRA IDs to process
        batch_size: Number of IDs to process in a single batch
        max_workers: Maximum number of parallel workers
        delay_between_batches: Delay in seconds between batches
        cache_file: Path to a JSON file to cache results (set to None to disable caching)
        max_retries: Maximum number of retries for failed requests
        retry_delay: Delay in seconds between retries
        show_progress: Whether to show a progress bar (works in Jupyter notebooks)
        
    Returns:
        Dictionary mapping SRA IDs to dictionaries containing 'bioproject_id' and 'geo_id'
    """
    start_time = time.time()
    logger.info(f"Processing {len(sra_ids)} SRA IDs")
    
    # If show_progress is True, wrap the convert_sra_ids function with tqdm
    if show_progress:
        try:
            # Initialize progress bar
            pbar = tqdm(total=len(sra_ids), desc="Processing SRA IDs")
            
            # Create a callback function to update the progress bar
            def update_progress(processed_count):
                pbar.update(processed_count - pbar.n)
            
            # Call the original function with the callback
            results = convert_sra_ids(
                sra_ids=sra_ids,
                batch_size=batch_size,
                max_workers=max_workers,
                delay_between_batches=delay_between_batches,
                cache_file=cache_file,
                max_retries=max_retries,
                retry_delay=retry_delay,
                progress_callback=update_progress
            )
            
            # Close the progress bar
            pbar.close()
        except Exception as e:
            logger.error(f"Error with progress bar: {str(e)}")
            # Fall back to no progress bar
            results = convert_sra_ids(
                sra_ids=sra_ids,
                batch_size=batch_size,
                max_workers=max_workers,
                delay_between_batches=delay_between_batches,
                cache_file=cache_file,
                max_retries=max_retries,
                retry_delay=retry_delay
            )
    else:
        # Use the original function without progress bar
        results = convert_sra_ids(
            sra_ids=sra_ids,
            batch_size=batch_size,
            max_workers=max_workers,
            delay_between_batches=delay_between_batches,
            cache_file=cache_file,
            max_retries=max_retries,
            retry_delay=retry_delay
        )
    
    elapsed_time = time.time() - start_time
    logger.info(f"Processing completed in {elapsed_time:.2f} seconds")
    
    # Print summary
    bioproject_count = sum(1 for result in results.values() if result['bioproject_id'])
    geo_count = sum(1 for result in results.values() if result['geo_id'])
    
    logger.info(f"Summary:")
    logger.info(f"  Total SRA IDs processed: {len(sra_ids)}")
    logger.info(f"  BioProject IDs found: {bioproject_count} ({bioproject_count/len(sra_ids)*100:.1f}%)")
    logger.info(f"  GEO/ArrayExpress IDs found: {geo_count} ({geo_count/len(sra_ids)*100:.1f}%)")
    
    return results

def process_sra_ids_from_file(
    file_path: str,
    output_file: Optional[str] = None,
    **kwargs
) -> Dict[str, Dict[str, str]]:
    """
    Process SRA IDs from a file and optionally save results to a TSV file.
    
    Args:
        file_path: Path to the file containing SRA IDs (one per line)
        output_file: Path to save the results as a TSV file (optional)
        **kwargs: Additional arguments to pass to process_sra_ids
        
    Returns:
        Dictionary mapping SRA IDs to dictionaries containing 'bioproject_id' and 'geo_id'
    """
    # Read SRA IDs from file
    logger.info(f"Reading SRA IDs from {file_path}")
    sra_ids = read_sra_ids_from_file(file_path)
    logger.info(f"Found {len(sra_ids)} SRA IDs in the input file")
    
    # Process SRA IDs
    results = process_sra_ids(sra_ids, **kwargs)
    
    # Save results to file if output_file is provided
    if output_file:
        logger.info(f"Saving results to {output_file}")
        results_to_tsv(results, output_file)
    
    return results

def results_to_dataframe(results: Dict[str, Dict[str, str]]) -> pd.DataFrame:
    """
    Convert results dictionary to a pandas DataFrame.
    
    Args:
        results: Dictionary mapping SRA IDs to dictionaries containing 'bioproject_id' and 'geo_id'
        
    Returns:
        pandas DataFrame with columns: SRA_ID, BioProject_ID, GEO_ArrayExpress_ID
    """
    data = []
    for sra_id, result in results.items():
        data.append({
            'SRA_ID': sra_id,
            'BioProject_ID': result.get('bioproject_id', ''),
            'GEO_ArrayExpress_ID': result.get('geo_id', '')
        })
    
    return pd.DataFrame(data)

def results_to_tsv(results: Dict[str, Dict[str, str]], output_file: str) -> None:
    """
    Save results to a TSV file.
    
    Args:
        results: Dictionary mapping SRA IDs to dictionaries containing 'bioproject_id' and 'geo_id'
        output_file: Path to the output TSV file
    """
    df = results_to_dataframe(results)
    df.to_csv(output_file, sep='\t', index=False)
    logger.info(f"Results saved to {output_file}")

# Example usage in a Jupyter notebook:
"""
from sra_batch_processor import process_sra_ids, process_sra_ids_from_file, results_to_dataframe

# Option 1: Process a list of SRA IDs
sra_ids = ['ERP149679', 'ERP144781', 'SRP324458', 'ERP127673']
results = process_sra_ids(
    sra_ids=sra_ids,
    batch_size=10,
    max_workers=2,
    delay_between_batches=3.0,
    cache_file="my_cache.json"
)

# Convert to DataFrame for easy viewing in notebook
df = results_to_dataframe(results)
display(df)

# Option 2: Process SRA IDs from a file
results = process_sra_ids_from_file(
    file_path="my_sra_ids.txt",
    output_file="results.tsv",
    batch_size=10,
    max_workers=2
)

# Display results
df = results_to_dataframe(results)
display(df)
""" 