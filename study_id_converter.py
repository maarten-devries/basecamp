import requests
import pandas as pd
import time
import re
from concurrent.futures import ThreadPoolExecutor
from typing import Dict, List, Tuple, Optional, Set
import xml.etree.ElementTree as ET
import os
import random

# NCBI API key - set this as an environment variable for better performance
# You can get an API key from: https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/
NCBI_API_KEY = os.environ.get('NCBI_API_KEY', '')

def convert_study_ids(study_ids: List[str], batch_size: int = 10, max_workers: int = 3, 
                     delay_between_batches: float = 1.0) -> Dict[str, str]:
    """
    Convert SRP/ERP IDs to their corresponding GSE/E-MTAB IDs.
    
    Args:
        study_ids: List of SRP or ERP IDs
        batch_size: Number of IDs to process in a single batch
        max_workers: Maximum number of parallel workers for API requests
        delay_between_batches: Delay in seconds between processing batches
        
    Returns:
        Dictionary mapping original IDs to their corresponding GSE/E-MTAB IDs
    """
    # Deduplicate study IDs
    unique_study_ids = list(set(study_ids))
    
    # Separate SRP and ERP IDs
    srp_ids = [sid for sid in unique_study_ids if sid.startswith('SRP')]
    erp_ids = [sid for sid in unique_study_ids if sid.startswith('ERP')]
    
    # Process in batches
    results = {}
    
    # Process SRP IDs (NCBI)
    if srp_ids:
        results.update(convert_srp_to_gse(srp_ids, batch_size, max_workers, delay_between_batches))
    
    # Process ERP IDs (EBI)
    if erp_ids:
        results.update(convert_erp_to_emtab(erp_ids, batch_size, max_workers, delay_between_batches))
    
    return results

def convert_srp_to_gse(srp_ids: List[str], batch_size: int = 10, max_workers: int = 3,
                      delay_between_batches: float = 1.0) -> Dict[str, str]:
    """
    Convert SRP IDs to GSE IDs using NCBI's Entrez API.
    
    Args:
        srp_ids: List of SRP IDs
        batch_size: Number of IDs to process in a single batch
        max_workers: Maximum number of parallel workers
        delay_between_batches: Delay in seconds between processing batches
        
    Returns:
        Dictionary mapping SRP IDs to GSE IDs
    """
    results = {}
    
    # Process in smaller batches to avoid overwhelming the API
    for i in range(0, len(srp_ids), batch_size):
        batch = srp_ids[i:i+batch_size]
        print(f"Processing SRP batch {i//batch_size + 1}/{(len(srp_ids) + batch_size - 1)//batch_size}: {batch}")
        
        try:
            batch_results = fetch_gse_for_srp_batch(batch)
            results.update(batch_results)
            
            # Count successful conversions in this batch
            successful = sum(1 for k, v in batch_results.items() if k != v)
            print(f"  Successfully converted {successful} out of {len(batch)} SRP IDs in this batch")
            
        except Exception as e:
            print(f"Error processing batch: {str(e)}")
            # If batch fails, try individual processing with longer delays
            print("  Trying individual processing with longer delays...")
            for srp_id in batch:
                try:
                    # Add a random delay to avoid hitting rate limits
                    time.sleep(random.uniform(0.5, 1.5))
                    individual_result = fetch_gse_for_srp_individual(srp_id)
                    results[srp_id] = individual_result
                    print(f"  {srp_id} -> {individual_result if individual_result != srp_id else '(not converted)'}")
                except Exception as e:
                    print(f"  Error processing {srp_id}: {str(e)}")
                    results[srp_id] = srp_id
        
        # Be nice to the API - add a random delay between batches
        if i + batch_size < len(srp_ids):
            delay = delay_between_batches + random.uniform(0, 1.0)
            print(f"  Waiting {delay:.2f} seconds before next batch...")
            time.sleep(delay)
    
    return results

def fetch_gse_for_srp_batch(srp_ids: List[str]) -> Dict[str, str]:
    """
    Fetch GSE IDs for a batch of SRP IDs using NCBI's Entrez API.
    
    Args:
        srp_ids: List of SRP IDs
        
    Returns:
        Dictionary mapping SRP IDs to GSE IDs
    """
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    results = {srp_id: srp_id for srp_id in srp_ids}  # Default to original IDs
    
    try:
        # Step 1: Search for SRP IDs in the SRA database
        search_url = f"{base_url}/esearch.fcgi"
        search_params = {
            "db": "sra",
            "term": " OR ".join(srp_ids),
            "retmode": "json",
            "retmax": 1000  # Increase if you have more than 1000 IDs in a batch
        }
        
        # Add API key if available
        if NCBI_API_KEY:
            search_params["api_key"] = NCBI_API_KEY
        
        search_response = requests.get(search_url, params=search_params)
        search_response.raise_for_status()
        search_data = search_response.json()
        
        # Get the UIDs from the search results
        uids = search_data.get("esearchresult", {}).get("idlist", [])
        
        if not uids:
            return results  # Return original IDs if not found
        
        # Step 2: Link to GEO database
        link_url = f"{base_url}/elink.fcgi"
        link_params = {
            "dbfrom": "sra",
            "db": "gds",  # GEO DataSets database
            "id": ",".join(uids),
            "retmode": "json"
        }
        
        # Add API key if available
        if NCBI_API_KEY:
            link_params["api_key"] = NCBI_API_KEY
        
        link_response = requests.get(link_url, params=link_params)
        link_response.raise_for_status()
        link_data = link_response.json()
        
        # Extract linked IDs
        linksets = link_data.get("linksets", [])
        
        # Create a mapping from SRA UIDs to SRP IDs
        # First, get summaries for the SRA UIDs
        summary_url = f"{base_url}/esummary.fcgi"
        summary_params = {
            "db": "sra",
            "id": ",".join(uids),
            "retmode": "json"
        }
        
        # Add API key if available
        if NCBI_API_KEY:
            summary_params["api_key"] = NCBI_API_KEY
        
        summary_response = requests.get(summary_url, params=summary_params)
        summary_response.raise_for_status()
        summary_data = summary_response.json()
        
        # Map SRA UIDs to SRP IDs
        uid_to_srp = {}
        for uid in uids:
            if uid in summary_data.get("result", {}):
                # Extract SRP ID from the summary
                expxml = summary_data["result"][uid].get("expxml", "")
                for srp_id in srp_ids:
                    if srp_id in expxml:
                        uid_to_srp[uid] = srp_id
                        break
        
        # Process linksets to get GEO IDs
        for linkset in linksets:
            sra_uid = linkset.get("ids", [None])[0]
            if not sra_uid or sra_uid not in uid_to_srp:
                continue
                
            srp_id = uid_to_srp[sra_uid]
            
            # Get linked GEO IDs
            linked_ids = []
            for linksetdb in linkset.get("linksetdbs", []):
                linked_ids.extend(linksetdb.get("links", []))
            
            if not linked_ids:
                continue
            
            # Fetch summary for the linked GEO datasets
            geo_summary_params = {
                "db": "gds",
                "id": ",".join(linked_ids),
                "retmode": "json"
            }
            
            # Add API key if available
            if NCBI_API_KEY:
                geo_summary_params["api_key"] = NCBI_API_KEY
            
            geo_summary_response = requests.get(summary_url, params=geo_summary_params)
            geo_summary_response.raise_for_status()
            geo_summary_data = geo_summary_response.json()
            
            # Extract GSE ID from summary
            for uid in linked_ids:
                if uid in geo_summary_data.get("result", {}):
                    accession = geo_summary_data["result"][uid].get("accession", "")
                    if accession.startswith("GSE"):
                        results[srp_id] = accession
                        break
                    
                    # If we couldn't find a GSE ID, try to extract it from the title or summary
                    title = geo_summary_data["result"][uid].get("title", "")
                    summary = geo_summary_data["result"][uid].get("summary", "")
                    
                    # Look for GSE pattern in title or summary
                    gse_match = re.search(r'(GSE\d+)', title + " " + summary)
                    if gse_match:
                        results[srp_id] = gse_match.group(1)
                        break
        
        return results
    
    except Exception as e:
        print(f"Error fetching GSE for batch: {str(e)}")
        raise  # Re-raise to handle in the calling function

def fetch_gse_for_srp_individual(srp_id: str) -> str:
    """
    Fetch GSE ID for a single SRP ID using NCBI's Entrez API.
    
    Args:
        srp_id: SRP ID
        
    Returns:
        GSE ID or original SRP ID if not found
    """
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    
    try:
        # Step 1: Search for the SRP ID in the SRA database
        search_url = f"{base_url}/esearch.fcgi"
        search_params = {
            "db": "sra",
            "term": srp_id,
            "retmode": "json"
        }
        
        # Add API key if available
        if NCBI_API_KEY:
            search_params["api_key"] = NCBI_API_KEY
        
        search_response = requests.get(search_url, params=search_params)
        search_response.raise_for_status()
        search_data = search_response.json()
        
        # Get the UIDs from the search results
        uids = search_data.get("esearchresult", {}).get("idlist", [])
        
        if not uids:
            return srp_id  # Return original ID if not found
        
        # Step 2: Link to GEO database
        link_url = f"{base_url}/elink.fcgi"
        link_params = {
            "dbfrom": "sra",
            "db": "gds",  # GEO DataSets database
            "id": uids[0],  # Just use the first UID
            "retmode": "json"
        }
        
        # Add API key if available
        if NCBI_API_KEY:
            link_params["api_key"] = NCBI_API_KEY
        
        link_response = requests.get(link_url, params=link_params)
        link_response.raise_for_status()
        link_data = link_response.json()
        
        # Extract linked IDs
        linksets = link_data.get("linksets", [])
        if not linksets or "linksetdbs" not in linksets[0]:
            return srp_id  # Return original ID if no links found
        
        linked_ids = []
        for linksetdb in linksets[0].get("linksetdbs", []):
            linked_ids.extend(linksetdb.get("links", []))
        
        if not linked_ids:
            return srp_id  # Return original ID if no links found
        
        # Step 3: Fetch summary for the linked GEO dataset
        summary_url = f"{base_url}/esummary.fcgi"
        summary_params = {
            "db": "gds",
            "id": linked_ids[0],  # Just use the first linked ID
            "retmode": "json"
        }
        
        # Add API key if available
        if NCBI_API_KEY:
            summary_params["api_key"] = NCBI_API_KEY
        
        summary_response = requests.get(summary_url, params=summary_params)
        summary_response.raise_for_status()
        summary_data = summary_response.json()
        
        # Extract GSE ID from summary
        if linked_ids[0] in summary_data.get("result", {}):
            accession = summary_data["result"][linked_ids[0]].get("accession", "")
            if accession.startswith("GSE"):
                return accession
            
            # If we couldn't find a GSE ID, try to extract it from the title or summary
            title = summary_data["result"][linked_ids[0]].get("title", "")
            summary_text = summary_data["result"][linked_ids[0]].get("summary", "")
            
            # Look for GSE pattern in title or summary
            gse_match = re.search(r'(GSE\d+)', title + " " + summary_text)
            if gse_match:
                return gse_match.group(1)
        
        # Step 4: Try a direct search in GEO as a fallback
        geo_search_params = {
            "db": "gds",
            "term": f"{srp_id}[Accession]",
            "retmode": "json"
        }
        
        # Add API key if available
        if NCBI_API_KEY:
            geo_search_params["api_key"] = NCBI_API_KEY
        
        geo_search_response = requests.get(search_url, params=geo_search_params)
        if geo_search_response.status_code == 200:
            geo_search_data = geo_search_response.json()
            geo_uids = geo_search_data.get("esearchresult", {}).get("idlist", [])
            
            if geo_uids:
                geo_summary_params = {
                    "db": "gds",
                    "id": geo_uids[0],
                    "retmode": "json"
                }
                
                # Add API key if available
                if NCBI_API_KEY:
                    geo_summary_params["api_key"] = NCBI_API_KEY
                
                geo_summary_response = requests.get(summary_url, params=geo_summary_params)
                if geo_summary_response.status_code == 200:
                    geo_summary_data = geo_summary_response.json()
                    
                    if geo_uids[0] in geo_summary_data.get("result", {}):
                        accession = geo_summary_data["result"][geo_uids[0]].get("accession", "")
                        if accession.startswith("GSE"):
                            return accession
        
        return srp_id  # Return original ID if GSE not found
    
    except Exception as e:
        print(f"Error fetching GSE for {srp_id}: {str(e)}")
        return srp_id  # Return original ID on error

def convert_erp_to_emtab(erp_ids: List[str], batch_size: int = 10, max_workers: int = 3,
                        delay_between_batches: float = 1.0) -> Dict[str, str]:
    """
    Convert ERP IDs to E-MTAB IDs using EBI's ArrayExpress API.
    
    Args:
        erp_ids: List of ERP IDs
        batch_size: Number of IDs to process in a single batch
        max_workers: Maximum number of parallel workers
        delay_between_batches: Delay in seconds between processing batches
        
    Returns:
        Dictionary mapping ERP IDs to E-MTAB IDs
    """
    results = {}
    
    # Process in batches to avoid overwhelming the API
    for i in range(0, len(erp_ids), batch_size):
        batch = erp_ids[i:i+batch_size]
        print(f"Processing ERP batch {i//batch_size + 1}/{(len(erp_ids) + batch_size - 1)//batch_size}: {batch}")
        
        # Use ThreadPoolExecutor for parallel processing
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            batch_results = list(executor.map(fetch_emtab_for_erp, batch))
        
        # Update results dictionary
        for erp_id, emtab_id in batch_results:
            results[erp_id] = emtab_id
        
        # Count successful conversions in this batch
        successful = sum(1 for erp_id, emtab_id in batch_results if erp_id != emtab_id)
        print(f"  Successfully converted {successful} out of {len(batch)} ERP IDs in this batch")
        
        # Be nice to the API - add a random delay between batches
        if i + batch_size < len(erp_ids):
            delay = delay_between_batches + random.uniform(0, 1.0)
            print(f"  Waiting {delay:.2f} seconds before next batch...")
            time.sleep(delay)
    
    return results

def fetch_emtab_for_erp(erp_id: str) -> Tuple[str, str]:
    """
    Fetch E-MTAB ID for a single ERP ID using EBI's ArrayExpress API.
    
    Args:
        erp_id: ERP ID
        
    Returns:
        Tuple of (ERP ID, E-MTAB ID or original ERP ID if not found)
    """
    # ArrayExpress API URL
    base_url = "https://www.ebi.ac.uk/arrayexpress/json/v3/experiments"
    
    try:
        # Search for the ERP ID in ArrayExpress
        search_url = f"{base_url}/{erp_id}"
        response = requests.get(search_url)
        
        # If direct lookup fails, try searching
        if response.status_code != 200:
            search_url = f"{base_url}?accession={erp_id}"
            response = requests.get(search_url)
        
        if response.status_code != 200:
            # Try another approach - search by secondary accession
            search_url = f"{base_url}?secondaryaccession={erp_id}"
            response = requests.get(search_url)
        
        if response.status_code == 200:
            data = response.json()
            experiments = data.get("experiments", {}).get("experiment", [])
            
            if experiments:
                # Look for E-MTAB ID in the accession field
                for exp in experiments:
                    accession = exp.get("accession", "")
                    if accession.startswith("E-MTAB-"):
                        return erp_id, accession
                    
                    # Check secondary accessions
                    sec_accessions = exp.get("secondaryaccession", "").split(" ")
                    for acc in sec_accessions:
                        if acc.startswith("E-MTAB-"):
                            return erp_id, acc
        
        # If we couldn't find an E-MTAB ID, try the BioStudies API
        biostudies_url = f"https://www.ebi.ac.uk/biostudies/api/v1/search?accession={erp_id}"
        biostudies_response = requests.get(biostudies_url)
        
        if biostudies_response.status_code == 200:
            biostudies_data = biostudies_response.json()
            hits = biostudies_data.get("hits", [])
            
            for hit in hits:
                accession = hit.get("accession", "")
                if accession.startswith("E-MTAB-"):
                    return erp_id, accession
        
        return erp_id, erp_id  # Return original ID if E-MTAB not found
    
    except Exception as e:
        print(f"Error fetching E-MTAB for {erp_id}: {str(e)}")
        return erp_id, erp_id  # Return original ID on error

def add_geo_emtab_ids_to_dataframe(df: pd.DataFrame, study_id_col: str = 'study_id', 
                                  batch_size: int = 10, max_workers: int = 3,
                                  delay_between_batches: float = 1.0) -> pd.DataFrame:
    """
    Add GSE/E-MTAB IDs to a dataframe based on SRP/ERP IDs.
    
    Args:
        df: Pandas DataFrame containing study IDs
        study_id_col: Column name containing the study IDs
        batch_size: Number of IDs to process in a single batch
        max_workers: Maximum number of parallel workers
        delay_between_batches: Delay in seconds between processing batches
        
    Returns:
        DataFrame with added 'geo_id' column
    """
    # Get unique study IDs
    unique_study_ids = df[study_id_col].unique().tolist()
    
    print(f"Converting {len(unique_study_ids)} unique study IDs...")
    
    # Convert IDs
    id_mapping = convert_study_ids(unique_study_ids, batch_size, max_workers, delay_between_batches)
    
    # Add new column to dataframe
    df['geo_id'] = df[study_id_col].map(id_mapping)
    
    # Count how many were successfully converted
    converted = sum(1 for k, v in id_mapping.items() if k != v)
    print(f"Successfully converted {converted} out of {len(unique_study_ids)} study IDs.")
    
    return df

# Example usage:
if __name__ == "__main__":
    # Example dataframe
    data = {'study_id': ['SRP559437', 'SRP270870', 'ERP149679']}
    df = pd.DataFrame(data)
    
    # Add GSE/E-MTAB IDs
    result_df = add_geo_emtab_ids_to_dataframe(df)
    print(result_df) 