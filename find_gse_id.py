import pandas as pd
import re
import requests
from typing import Optional, Dict, Any, Union

def find_gse_id(entrez_id=None, srx_accession=None) -> Optional[str]:
    """
    Find the GSE ID for a given entrez_id or SRX accession.
    
    This is a simple function that uses a direct approach to find GSE IDs.
    It first tries to use the ENA API for ERX accessions, then falls back to
    a simple regex search in the NCBI Entrez API response.
    
    Args:
        entrez_id: The entrez_id to search for
        srx_accession: The SRX accession to search for
        
    Returns:
        The GSE ID if found, None otherwise
    """
    # Try using SRX accession first
    if srx_accession:
        # For ERX accessions, try the ENA API
        if srx_accession.startswith('ERX'):
            try:
                # ENA API URL
                ena_url = "https://www.ebi.ac.uk/ena/browser/api/xml"
                url = f"{ena_url}/{srx_accession}"
                
                response = requests.get(url)
                response.raise_for_status()
                
                # Look for GSE ID in the XML content
                content = response.text
                gse_match = re.search(r'(GSE\d+)', content)
                
                if gse_match:
                    return gse_match.group(1)
            except Exception as e:
                print(f"Error querying ENA API: {e}")
        
        # Try NCBI Entrez API
        try:
            # Base URL for NCBI Entrez API
            base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
            
            # Search for the SRX accession in the SRA database
            search_url = f"{base_url}esearch.fcgi"
            search_params = {
                "db": "sra",
                "term": srx_accession,
                "retmax": 1
            }
            
            search_response = requests.get(search_url, params=search_params)
            search_response.raise_for_status()
            
            # Get the ID from the search result
            search_content = search_response.text
            id_match = re.search(r'<Id>(\d+)</Id>', search_content)
            
            if id_match:
                sra_id = id_match.group(1)
                
                # Use the SRA ID to get the full record
                fetch_url = f"{base_url}efetch.fcgi"
                fetch_params = {
                    "db": "sra",
                    "id": sra_id,
                    "retmode": "xml"
                }
                
                fetch_response = requests.get(fetch_url, params=fetch_params)
                fetch_response.raise_for_status()
                
                # Look for GSE ID in the XML content
                fetch_content = fetch_response.text
                gse_match = re.search(r'(GSE\d+)', fetch_content)
                
                if gse_match:
                    return gse_match.group(1)
        except Exception as e:
            print(f"Error querying NCBI API: {e}")
    
    # Try using entrez_id
    if entrez_id:
        try:
            # Base URL for NCBI Entrez API
            base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
            
            # Search for the entrez_id in the SRA database
            search_url = f"{base_url}esearch.fcgi"
            search_params = {
                "db": "sra",
                "term": str(entrez_id),
                "retmax": 1
            }
            
            search_response = requests.get(search_url, params=search_params)
            search_response.raise_for_status()
            
            # Get the ID from the search result
            search_content = search_response.text
            id_match = re.search(r'<Id>(\d+)</Id>', search_content)
            
            if id_match:
                sra_id = id_match.group(1)
                
                # Use the SRA ID to get the full record
                fetch_url = f"{base_url}efetch.fcgi"
                fetch_params = {
                    "db": "sra",
                    "id": sra_id,
                    "retmode": "xml"
                }
                
                fetch_response = requests.get(fetch_url, params=fetch_params)
                fetch_response.raise_for_status()
                
                # Look for GSE ID in the XML content
                fetch_content = fetch_response.text
                gse_match = re.search(r'(GSE\d+)', fetch_content)
                
                if gse_match:
                    return gse_match.group(1)
        except Exception as e:
            print(f"Error querying NCBI API: {e}")
    
    return None

def add_gse_ids_to_df(df, entrez_id_col='entrez_id', srx_col='srx_accession'):
    """
    Add GSE IDs to a dataframe based on entrez_ids and SRX accessions.
    
    Args:
        df: Pandas DataFrame containing entrez_ids and SRX accessions
        entrez_id_col: Name of the column containing entrez_ids
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
            
        entrez_id = row.get(entrez_id_col)
        srx_accession = row.get(srx_col)
        
        # Try to get the GSE ID
        gse_id = find_gse_id(entrez_id=entrez_id, srx_accession=srx_accession)
        
        # Update the dataframe
        result_df.loc[idx, 'gse_id'] = gse_id
        
        # Print progress
        if (i + 1) % 10 == 0 or i == len(result_df) - 1:
            print(f"Processed {i + 1}/{len(result_df)} entries")
    
    return result_df

# Example usage
if __name__ == "__main__":
    # Test with an example
    test_entrez_id = 29110018
    test_srx = "ERX11148735"
    
    print(f"Testing with entrez_id: {test_entrez_id}")
    gse_id = find_gse_id(entrez_id=test_entrez_id)
    print(f"GSE ID: {gse_id}")
    
    print(f"Testing with SRX accession: {test_srx}")
    gse_id = find_gse_id(srx_accession=test_srx)
    print(f"GSE ID: {gse_id}")
    
    # Test with a dataframe
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