#!/usr/bin/env python3
import argparse
import pandas as pd
from add_gse_ids import add_gse_ids_to_df, add_gse_ids_to_df_parallel

def main():
    parser = argparse.ArgumentParser(description='Add GSE IDs to a CSV file based on entrez_ids and SRX accessions.')
    parser.add_argument('input_file', help='Input CSV file')
    parser.add_argument('--output_file', help='Output CSV file (default: input_file with _with_gse_ids suffix)')
    parser.add_argument('--entrez_id_col', default='entrez_id', help='Name of the column containing entrez_ids')
    parser.add_argument('--srx_col', default='srx_accession', help='Name of the column containing SRX accessions')
    parser.add_argument('--batch_size', type=int, default=10, help='Number of rows to process before saving intermediate results')
    parser.add_argument('--parallel', action='store_true', help='Use parallel processing')
    parser.add_argument('--n_jobs', type=int, default=4, help='Number of parallel jobs')
    
    args = parser.parse_args()
    
    # Set default output file if not provided
    if args.output_file is None:
        args.output_file = args.input_file.replace('.csv', '_with_gse_ids.csv')
        if args.output_file == args.input_file:
            args.output_file = args.input_file + '_with_gse_ids.csv'
    
    # Load the input file
    print(f"Loading {args.input_file}...")
    df = pd.read_csv(args.input_file)
    print(f"Loaded {len(df)} rows.")
    
    # Check if the required columns exist
    if args.entrez_id_col not in df.columns:
        print(f"Error: Column '{args.entrez_id_col}' not found in the input file.")
        return
    
    if args.srx_col not in df.columns:
        print(f"Error: Column '{args.srx_col}' not found in the input file.")
        return
    
    # Process the dataframe
    print(f"Processing dataframe to add GSE IDs...")
    if args.parallel:
        print(f"Using parallel processing with {args.n_jobs} jobs.")
        result_df = add_gse_ids_to_df_parallel(df, args.entrez_id_col, args.srx_col, args.n_jobs)
    else:
        print(f"Using sequential processing with batch size {args.batch_size}.")
        result_df = add_gse_ids_to_df(df, args.entrez_id_col, args.srx_col, args.batch_size)
    
    # Save the result
    print(f"Saving result to {args.output_file}...")
    result_df.to_csv(args.output_file, index=False)
    print(f"Done! Result saved to {args.output_file}.")
    
    # Print summary
    total_rows = len(result_df)
    found_gse_ids = result_df['gse_id'].notna().sum()
    print(f"Summary: Found GSE IDs for {found_gse_ids} out of {total_rows} rows ({found_gse_ids/total_rows*100:.2f}%).")

if __name__ == "__main__":
    main() 