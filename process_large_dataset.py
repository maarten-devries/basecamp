#!/usr/bin/env python3
"""
Script to demonstrate processing a large list of SRA IDs efficiently.
This script shows how to use the enhanced convert_sra_ids function with caching
and batching to handle large datasets without hitting API limits.
"""

import os
import sys
import logging
import time
import argparse
from typing import List
from sra_id_converter import convert_sra_ids

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

def write_results_to_file(results, output_file: str):
    """
    Write conversion results to a tab-separated file.
    
    Args:
        results: Dictionary mapping SRA IDs to their conversion results
        output_file: Path to the output file
    """
    with open(output_file, 'w') as f:
        # Write header
        f.write("SRA_ID\tBioProject_ID\tGEO_ArrayExpress_ID\n")
        
        # Write results
        for sra_id, result in results.items():
            f.write(f"{sra_id}\t{result['bioproject_id']}\t{result['geo_id']}\n")

def main():
    parser = argparse.ArgumentParser(description="Convert a large list of SRA IDs to BioProject and GEO/ArrayExpress IDs")
    parser.add_argument("input_file", help="Path to file containing SRA IDs (one per line)")
    parser.add_argument("--output", "-o", default="sra_conversion_results.tsv", help="Path to output TSV file")
    parser.add_argument("--cache", "-c", default="sra_cache.json", help="Path to cache file")
    parser.add_argument("--batch-size", "-b", type=int, default=10, help="Number of IDs to process in a single batch")
    parser.add_argument("--max-workers", "-w", type=int, default=2, help="Maximum number of parallel workers")
    parser.add_argument("--delay", "-d", type=float, default=3.0, help="Delay in seconds between batches")
    parser.add_argument("--max-retries", "-r", type=int, default=3, help="Maximum number of retries for failed requests")
    
    args = parser.parse_args()
    
    # Read SRA IDs from file
    logger.info(f"Reading SRA IDs from {args.input_file}")
    try:
        sra_ids = read_sra_ids_from_file(args.input_file)
    except Exception as e:
        logger.error(f"Error reading input file: {str(e)}")
        sys.exit(1)
    
    logger.info(f"Found {len(sra_ids)} SRA IDs in the input file")
    
    # Process SRA IDs
    start_time = time.time()
    logger.info("Starting SRA ID conversion")
    
    try:
        results = convert_sra_ids(
            sra_ids=sra_ids,
            batch_size=args.batch_size,
            max_workers=args.max_workers,
            delay_between_batches=args.delay,
            cache_file=args.cache,
            max_retries=args.max_retries
        )
    except Exception as e:
        logger.error(f"Error during conversion: {str(e)}")
        sys.exit(1)
    
    elapsed_time = time.time() - start_time
    logger.info(f"Conversion completed in {elapsed_time:.2f} seconds")
    
    # Write results to file
    logger.info(f"Writing results to {args.output}")
    try:
        write_results_to_file(results, args.output)
    except Exception as e:
        logger.error(f"Error writing output file: {str(e)}")
        sys.exit(1)
    
    # Print summary
    bioproject_count = sum(1 for result in results.values() if result['bioproject_id'])
    geo_count = sum(1 for result in results.values() if result['geo_id'])
    
    logger.info(f"Summary:")
    logger.info(f"  Total SRA IDs processed: {len(sra_ids)}")
    logger.info(f"  BioProject IDs found: {bioproject_count} ({bioproject_count/len(sra_ids)*100:.1f}%)")
    logger.info(f"  GEO/ArrayExpress IDs found: {geo_count} ({geo_count/len(sra_ids)*100:.1f}%)")
    logger.info(f"Results written to {args.output}")

if __name__ == "__main__":
    main() 