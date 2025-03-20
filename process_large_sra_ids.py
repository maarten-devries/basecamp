#!/usr/bin/env python3
"""
Example script demonstrating how to use the sra_batch_processor module
to process a large list of SRA IDs.
"""

import argparse
import logging
import time
from sra_batch_processor import process_sra_ids_from_file, results_to_dataframe

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("sra_conversion.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

def main():
    parser = argparse.ArgumentParser(description="Process a large list of SRA IDs")
    parser.add_argument("input_file", help="Path to file containing SRA IDs (one per line)")
    parser.add_argument("--output", "-o", default="sra_conversion_results.tsv", help="Path to output TSV file")
    parser.add_argument("--cache", "-c", default="sra_cache.json", help="Path to cache file")
    parser.add_argument("--batch-size", "-b", type=int, default=20, help="Number of IDs to process in a single batch")
    parser.add_argument("--max-workers", "-w", type=int, default=3, help="Maximum number of parallel workers")
    parser.add_argument("--delay", "-d", type=float, default=3.0, help="Delay in seconds between batches")
    parser.add_argument("--max-retries", "-r", type=int, default=3, help="Maximum number of retries for failed requests")
    parser.add_argument("--no-progress", action="store_true", help="Disable progress bar")
    
    args = parser.parse_args()
    
    start_time = time.time()
    logger.info(f"Starting SRA ID conversion from file: {args.input_file}")
    
    # Process SRA IDs from file
    results = process_sra_ids_from_file(
        file_path=args.input_file,
        output_file=args.output,
        batch_size=args.batch_size,
        max_workers=args.max_workers,
        delay_between_batches=args.delay,
        cache_file=args.cache,
        max_retries=args.max_retries,
        show_progress=not args.no_progress
    )
    
    elapsed_time = time.time() - start_time
    logger.info(f"Total processing time: {elapsed_time:.2f} seconds")
    
    # Print summary
    df = results_to_dataframe(results)
    logger.info(f"Results summary:")
    logger.info(f"  Total SRA IDs: {len(df)}")
    logger.info(f"  SRA IDs with BioProject ID: {df['BioProject_ID'].notna().sum()} ({df['BioProject_ID'].notna().mean()*100:.1f}%)")
    logger.info(f"  SRA IDs with GEO/ArrayExpress ID: {df['GEO_ArrayExpress_ID'].notna().sum()} ({df['GEO_ArrayExpress_ID'].notna().mean()*100:.1f}%)")
    logger.info(f"Results saved to: {args.output}")

if __name__ == "__main__":
    main() 