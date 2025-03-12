    main() 
#!/usr/bin/env python3
"""
A simple script to find GSE IDs for SRX/ERX accessions.
"""

import pandas as pd
import re
import sys
import argparse
import requests
from typing import Optional

def get_gse_id_from_ena(srx_accession: str) -> Optional[str]:
    """
    Find the GSE ID for a given SRX/ERX accession by querying the ENA API.
    
    Args:
        srx_accession: The SRX/ERX accession to search for
        
    Returns:
        The GSE ID if found, None otherwise
    """
    # ENA API URL
    ena_url = "https://www.ebi.ac.uk/ena/browser/api/xml"
    
    try:
        # Make the API request
        url = f"{ena_url}/{srx_accession}"
        print(f"Making request to {url}")
        response = requests.get(url)
        response.raise_for_status()
        
        # Print the response for debugging
        print(f"Response status code: {response.status_code}")
        
        # Look for GSE ID in the XML content
        content = response.text
        gse_match = re.search(r'(GSE\d+)', content)
        
        if gse_match:
            gse_id = gse_match.group(1)
            print(f"Found GSE ID: {gse_id}")
            return gse_id
        
        print(f"No GSE ID found for accession: {srx_accession}")
        return None
        
    except requests.exceptions.RequestException as e:
        print(f"Error querying ENA API: {e}")
        return None

def main():
    parser = argparse.ArgumentParser(description='Find GSE IDs for SRX/ERX accessions.')
    parser.add_argument('srx_accession', help='SRX accession to search for')
    
    args = parser.parse_args()
    
    # Process a single SRX accession
    gse_id = get_gse_id_from_ena(args.srx_accession)
    print(f"GSE ID for {args.srx_accession}: {gse_id}")

if __name__ == "__main__":
 
A simple script to find GSE IDs for SRX/ERX accessions.
"""

import pandas as pd
import re
import sys
import argparse
import requests
from typing import Optional, Dict, List, Union

def get_gse_id_from_ena(srx_accession: str) -> Optional[str]:
    """
    Find the GSE ID for a given SRX/ERX accession by querying the ENA API.
    
    Args:
        srx_accession: The SRX/ERX accession to search for
        
    Returns:
        The GSE ID if found, None otherwise
    """
    # ENA API URL
    ena_url = "https://www.ebi.ac.uk/ena/browser/api/xml"
    
    try:
        # Make the API request
        url = f"{ena_url}/{srx_accession}"
        print(f"Making request to {url}")
        response = requests.get(url)
        response.raise_for_status()
        
        # Print the response for debugging
        print(f"Response status code: {response.status_code}")
        
        # Look for GSE ID in the XML content
        content = response.text
        gse_match = re.search(r'(GSE\d+)', content)
        
        if gse_match:
            gse_id = gse_match.group(1)
            print(f"Found GSE ID: {gse_id}")
            return gse_id
        
        print(f"No GSE ID found for accession: {srx_accession}")
        return None
        
    except requests.exceptions.RequestException as e:
        print(f"Error querying ENA API: {e}")
        return None

def add_gse_ids_to_df(df: pd.DataFrame, srx_col: str = 'srx_accession') -> pd.DataFrame:
    """
    Add GSE IDs to a dataframe based on SRX accessions.
    
    Args:
        df: Pandas DataFrame containing SRX accessions
        srx_col: Name of the column containing SRX accessions
        
    Returns:
        DataFrame with an additional 'gse_id' column
    """
    # Create a copy of the dataframe to avoid modifying the original
    result_df = df.copy()
    
    # Add a new column for GSE IDs
    if 'gse_id' not in result_df.columns:
        result_df['gse_id'] = None
    
    # Process each row
    for i, (idx, row) in enumerate(result_df.iterrows()):
        # Skip if we already have a GSE ID for this row
        if pd.notna(result_df.loc[idx, 'gse_id']):
            continue
            
        srx_accession = row[srx_col]
        
        # Try to get the GSE ID
        gse_id = get_gse_id_from_ena(srx_accession)
        
        # Update the dataframe
        result_df.loc[idx, 'gse_id'] = gse_id
        
        # Print progress
        print(f"Processed {i + 1}/{len(result_df)} entries")
    
    return result_df

def main():
    parser = argparse.ArgumentParser(description='Find GSE IDs for SRX/ERX accessions.')
    parser.add_argument('--input_file', help='Input CSV file containing SRX accessions')
    parser.add_argument('--output_file', help='Output CSV file')
    parser.add_argument('--srx_col', default='srx_accession', help='Name of the column containing SRX accessions')
    parser.add_argument('srx_accession', nargs='?', help='SRX accession to search for')
    
    args = parser.parse_args()
    
    if args.input_file:
        # Process a CSV file
        df = pd.read_csv(args.input_file)
        result_df = add_gse_ids_to_df(df, args.srx_col)
        
        # Save the result
        output_file = args.output_file or args.input_file.replace('.csv', '_with_gse_ids.csv')
        result_df.to_csv(output_file, index=False)
        print(f"Result saved to {output_file}")
    elif args.srx_accession:
        # Process a single SRX accession
        gse_id = get_gse_id_from_ena(args.srx_accession)
        print(f"GSE ID for {args.srx_accession}: {gse_id}")
    else:
        parser.print_help()

if __name__ == "__main__":
 