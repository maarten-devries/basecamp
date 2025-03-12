import pandas as pd
import time
from get_gse_id import get_gse_id

def add_gse_ids_to_df(df, entrez_id_col='entrez_id', srx_col='srx_accession', batch_size=10):
    """
    Add GSE IDs to a dataframe based on entrez_ids and SRX accessions.
    
    Args:
        df: Pandas DataFrame containing entrez_ids and SRX accessions
        entrez_id_col: Name of the column containing entrez_ids
        srx_col: Name of the column containing SRX accessions
        batch_size: Number of rows to process before saving intermediate results
        
    Returns:
        DataFrame with an additional 'gse_id' column
    """
    # Create a copy of the dataframe to avoid modifying the original
    result_df = df.copy()
    
    # Add a new column for GSE IDs
    if 'gse_id' not in result_df.columns:
        result_df['gse_id'] = None
    
    # Process in batches to allow for intermediate saving
    total_rows = len(result_df)
    
    for i, (idx, row) in enumerate(result_df.iterrows()):
        # Skip if we already have a GSE ID for this row
        if pd.notna(result_df.loc[idx, 'gse_id']):
            continue
            
        entrez_id = row[entrez_id_col]
        srx_accession = row[srx_col]
        
        # Use the combined function to get the GSE ID
        gse_id = get_gse_id(entrez_id=entrez_id, srx_accession=srx_accession)
        
        # Update the dataframe
        result_df.loc[idx, 'gse_id'] = gse_id
        
        # Print progress
        if (i + 1) % 5 == 0 or i == total_rows - 1:
            print(f"Processed {i + 1}/{total_rows} entries")
        
        # Save intermediate results
        if (i + 1) % batch_size == 0:
            print(f"Saving intermediate results after processing {i + 1} entries")
            result_df.to_csv(f"df_with_gse_ids_intermediate_{i+1}.csv", index=False)
            
        # Add a small delay to avoid hitting NCBI rate limits
        time.sleep(0.3)
    
    # Save final results
    result_df.to_csv("df_with_gse_ids.csv", index=False)
    
    return result_df

def add_gse_ids_to_df_parallel(df, entrez_id_col='entrez_id', srx_col='srx_accession', n_jobs=4):
    """
    Add GSE IDs to a dataframe using parallel processing.
    
    Args:
        df: Pandas DataFrame containing entrez_ids and SRX accessions
        entrez_id_col: Name of the column containing entrez_ids
        srx_col: Name of the column containing SRX accessions
        n_jobs: Number of parallel jobs
        
    Returns:
        DataFrame with an additional 'gse_id' column
    """
    try:
        from joblib import Parallel, delayed
    except ImportError:
        print("joblib is required for parallel processing. Install it with 'pip install joblib'")
        return add_gse_ids_to_df(df, entrez_id_col, srx_col)
    
    # Create a copy of the dataframe to avoid modifying the original
    result_df = df.copy()
    
    # Add a new column for GSE IDs
    if 'gse_id' not in result_df.columns:
        result_df['gse_id'] = None
    
    def process_row(row):
        entrez_id = row[entrez_id_col]
        srx_accession = row[srx_col]
        
        # Use the combined function to get the GSE ID
        gse_id = get_gse_id(entrez_id=entrez_id, srx_accession=srx_accession)
        
        return gse_id
    
    # Process rows in parallel
    gse_ids = Parallel(n_jobs=n_jobs)(delayed(process_row)(row) for _, row in df.iterrows())
    
    # Update the dataframe
    result_df['gse_id'] = gse_ids
    
    # Save results
    result_df.to_csv("df_with_gse_ids.csv", index=False)
    
    return result_df

if __name__ == "__main__":
    # Example usage
    # Load your dataframe
    # df = pd.read_csv("your_data.csv")
    
    # For testing with a small example
    df = pd.DataFrame({
        'entrez_id': [29110018, 29110027],
        'srx_accession': ['ERX11148735', 'ERX11148744']
    })
    
    print("Original dataframe:")
    print(df)
    
    # Add GSE IDs
    result_df = add_gse_ids_to_df(df)
    
    print("\nDataframe with GSE IDs:")
    print(result_df) 