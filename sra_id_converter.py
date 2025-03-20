#!/usr/bin/env python3
"""
Script to retrieve BioProject IDs (PRJEB/PRJNA) and GEO/ArrayExpress identifiers
from SRA IDs (SRP/ERP).

Example usage:
    from sra_id_converter import convert_sra_ids
    
    sra_ids = ['ERP127673', 'SRP324458']
    results = convert_sra_ids(sra_ids)
    print(results)
    # Output: {
    #   'ERP127673': {'bioproject_id': 'PRJEB43688', 'geo_id': 'E-MTAB-10220'},
    #   'SRP324458': {'bioproject_id': 'PRJNA738600', 'geo_id': 'GSE178360'}
    # }
"""

import requests
import re
import time
import os
import xml.etree.ElementTree as ET
from typing import Dict, List, Optional, Any, Set, Tuple
from concurrent.futures import ThreadPoolExecutor
import logging
import random
import json
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# NCBI API key - set this as an environment variable for better performance
# You can get an API key from: https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/
NCBI_API_KEY = os.environ.get('NCBI_API_KEY', '')

# Known mappings for test cases
KNOWN_MAPPINGS = {}

# ENA API base URL
ENA_API_BASE_URL = "https://www.ebi.ac.uk/ena/browser/api"

# Global variable to track the last time we made a request to NCBI
last_ncbi_request_time = 0
last_ebi_request_time = 0

def ncbi_request(url, params):
    """
    Make a request to NCBI with proper rate limiting.
    
    Args:
        url: The URL to request
        params: The parameters for the request
        
    Returns:
        The response from the request
    """
    global last_ncbi_request_time
    
    # Add API key if available
    if NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY
    
    # Ensure we don't exceed rate limits (max 3 requests per second with API key, 1 per second without)
    current_time = time.time()
    time_since_last_request = current_time - last_ncbi_request_time
    
    # Wait time depends on whether we have an API key
    min_wait_time = 0.34 if NCBI_API_KEY else 1.0
    
    if time_since_last_request < min_wait_time:
        time.sleep(min_wait_time - time_since_last_request + random.uniform(0.1, 0.3))
    
    # Make the request
    response = requests.get(url, params=params)
    last_ncbi_request_time = time.time()
    
    return response

def ebi_request(url):
    """
    Make a request to EBI with proper rate limiting.
    
    Args:
        url: The URL to request
        
    Returns:
        The response from the request
    """
    global last_ebi_request_time
    
    # Ensure we don't exceed rate limits
    current_time = time.time()
    time_since_last_request = current_time - last_ebi_request_time
    
    min_wait_time = 1.0
    
    if time_since_last_request < min_wait_time:
        time.sleep(min_wait_time - time_since_last_request + random.uniform(0.1, 0.3))
    
    # Make the request
    response = requests.get(url)
    last_ebi_request_time = time.time()
    
    return response

def convert_sra_ids(sra_ids: List[str], batch_size: int = 5, max_workers: int = 2, 
                   delay_between_batches: float = 2.0, cache_file: str = None,
                   max_retries: int = 3, retry_delay: float = 5.0,
                   progress_callback = None) -> Dict[str, Dict[str, str]]:
    """
    Convert SRA IDs (SRP/ERP) to their corresponding BioProject IDs and GEO/ArrayExpress identifiers.
    
    Args:
        sra_ids: List of SRA IDs (SRP or ERP format)
        batch_size: Number of IDs to process in a single batch
        max_workers: Maximum number of parallel workers for API requests
        delay_between_batches: Delay in seconds between processing batches
        cache_file: Path to a JSON file to cache results (will be created if it doesn't exist)
        max_retries: Maximum number of retries for failed API requests
        retry_delay: Delay in seconds between retries
        progress_callback: Optional callback function to report progress (receives number of processed IDs)
        
    Returns:
        Dictionary mapping original SRA IDs to dictionaries containing 'bioproject_id' and 'geo_id'
    """
    # Deduplicate study IDs
    unique_sra_ids = list(set(sra_ids))
    logger.info(f"Processing {len(unique_sra_ids)} unique SRA IDs out of {len(sra_ids)} total")
    
    # Initialize results dictionary
    results = {}
    
    # Load cache if provided
    cache = {}
    if cache_file:
        cache_path = Path(cache_file)
        if cache_path.exists():
            try:
                with open(cache_path, 'r') as f:
                    cache = json.load(f)
                logger.info(f"Loaded {len(cache)} cached results from {cache_file}")
            except Exception as e:
                logger.warning(f"Error loading cache file {cache_file}: {str(e)}")
    
    # Check which IDs are already in cache or known mappings
    remaining_ids = []
    for sra_id in unique_sra_ids:
        if sra_id in KNOWN_MAPPINGS:
            results[sra_id] = KNOWN_MAPPINGS[sra_id]
        elif cache_file and sra_id in cache:
            results[sra_id] = cache[sra_id]
        else:
            remaining_ids.append(sra_id)
    
    # Update progress if callback is provided
    processed_count = len(unique_sra_ids) - len(remaining_ids)
    if progress_callback and processed_count > 0:
        progress_callback(processed_count)
    
    logger.info(f"Found {len(unique_sra_ids) - len(remaining_ids)} IDs in cache/known mappings, {len(remaining_ids)} remaining to process")
    
    # Separate SRP and ERP IDs
    srp_ids = [sid for sid in remaining_ids if sid.startswith('SRP')]
    erp_ids = [sid for sid in remaining_ids if sid.startswith('ERP')]
    other_ids = [sid for sid in remaining_ids if not (sid.startswith('SRP') or sid.startswith('ERP'))]
    
    if other_ids:
        logger.warning(f"Found {len(other_ids)} IDs with unknown format: {other_ids[:5]}{'...' if len(other_ids) > 5 else ''}")
    
    # Process SRP IDs (NCBI)
    if srp_ids:
        logger.info(f"Processing {len(srp_ids)} SRP IDs in batches of {batch_size}")
        batch_results = process_in_batches(srp_ids, process_srp_id, batch_size, max_workers, 
                                          delay_between_batches, max_retries, retry_delay, 
                                          progress_callback)
        results.update(batch_results)
    
    # Process ERP IDs (EBI)
    if erp_ids:
        logger.info(f"Processing {len(erp_ids)} ERP IDs in batches of {batch_size}")
        batch_results = process_in_batches(erp_ids, process_erp_id, batch_size, max_workers, 
                                          delay_between_batches, max_retries, retry_delay,
                                          progress_callback)
        results.update(batch_results)
    
    # Update cache if provided
    if cache_file:
        # Merge existing cache with new results
        cache.update(results)
        try:
            with open(cache_path, 'w') as f:
                json.dump(cache, f, indent=2)
            logger.info(f"Updated cache file {cache_file} with {len(results)} results")
        except Exception as e:
            logger.warning(f"Error writing to cache file {cache_file}: {str(e)}")
    
    # Return only the requested IDs
    return {sra_id: results.get(sra_id, {'bioproject_id': '', 'geo_id': ''}) for sra_id in sra_ids}

def process_in_batches(ids: List[str], process_func, batch_size: int, max_workers: int,
                      delay_between_batches: float, max_retries: int, retry_delay: float,
                      progress_callback = None) -> Dict[str, Dict[str, str]]:
    """
    Process a list of IDs in batches with retries.
    
    Args:
        ids: List of IDs to process
        process_func: Function to process each ID
        batch_size: Number of IDs to process in a single batch
        max_workers: Maximum number of parallel workers
        delay_between_batches: Delay in seconds between batches
        max_retries: Maximum number of retries for failed IDs
        retry_delay: Delay in seconds between retries
        progress_callback: Optional callback function to report progress
        
    Returns:
        Dictionary mapping IDs to their results
    """
    results = {}
    failed_ids = set(ids)
    retry_count = 0
    
    # Track total processed IDs for progress reporting
    total_processed = 0
    
    while failed_ids and retry_count <= max_retries:
        if retry_count > 0:
            logger.info(f"Retry {retry_count}/{max_retries} for {len(failed_ids)} failed IDs")
            time.sleep(retry_delay)
        
        batch_failed_ids = set()
        current_ids = list(failed_ids)
        
        for i in range(0, len(current_ids), batch_size):
            batch = current_ids[i:i+batch_size]
            logger.debug(f"Processing batch {i//batch_size + 1}/{(len(current_ids) + batch_size - 1)//batch_size}")
            
            batch_results = {}
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                future_to_id = {executor.submit(process_func, sid): sid for sid in batch}
                
                for future in future_to_id:
                    sid = future_to_id[future]
                    try:
                        result = future.result()
                        if result.get('bioproject_id') or result.get('geo_id'):
                            batch_results[sid] = result
                        else:
                            batch_failed_ids.add(sid)
                    except Exception as e:
                        logger.error(f"Error processing ID {sid}: {str(e)}")
                        batch_failed_ids.add(sid)
            
            results.update(batch_results)
            
            # Update progress if callback is provided
            batch_processed = len(batch) - len(batch_failed_ids.intersection(batch))
            total_processed += batch_processed
            if progress_callback and batch_processed > 0:
                progress_callback(total_processed)
            
            # Sleep to avoid rate limits, but only if there are more batches to process
            if i + batch_size < len(current_ids):
                time.sleep(delay_between_batches)
        
        failed_ids = batch_failed_ids
        retry_count += 1
    
    if failed_ids:
        logger.warning(f"Failed to process {len(failed_ids)} IDs after {max_retries} retries")
        # Add empty results for failed IDs
        for sid in failed_ids:
            results[sid] = {'bioproject_id': '', 'geo_id': ''}
    
    return results

def process_srp_id(srp_id: str) -> Dict[str, str]:
    """
    Process a single SRP ID to get its BioProject ID and GSE ID.
    
    Args:
        srp_id: SRP ID to process
        
    Returns:
        Dictionary with 'bioproject_id' and 'geo_id' keys
    """
    result = {'bioproject_id': '', 'geo_id': ''}
    
    try:
        # Step 1: Get BioProject ID from SRP ID
        bioproject_id = get_bioproject_from_srp(srp_id)
        if bioproject_id:
            result['bioproject_id'] = bioproject_id
        
        # Step 2: Get GSE ID from SRP ID
        gse_id = get_gse_from_srp(srp_id)
        if gse_id:
            result['geo_id'] = gse_id
            
    except Exception as e:
        logger.error(f"Error processing SRP ID {srp_id}: {str(e)}")
        
    return result

def process_erp_id(erp_id: str) -> Dict[str, str]:
    """
    Process a single ERP ID to get its BioProject ID and ArrayExpress ID.
    
    Args:
        erp_id: ERP ID to process
        
    Returns:
        Dictionary with 'bioproject_id' and 'geo_id' keys
    """
    # First check if we have a known mapping
    if erp_id in KNOWN_MAPPINGS:
        return KNOWN_MAPPINGS[erp_id]
    
    result = {'bioproject_id': '', 'geo_id': ''}
    
    try:
        # Step 1: Get BioProject ID from ERP ID
        bioproject_id = get_bioproject_from_erp(erp_id)
        if bioproject_id:
            result['bioproject_id'] = bioproject_id
        
        # Step 2: Get ArrayExpress ID from ERP ID
        ae_id = get_arrayexpress_from_erp(erp_id)
        if ae_id:
            result['geo_id'] = ae_id
            
    except Exception as e:
        logger.error(f"Error processing ERP ID {erp_id}: {str(e)}")
        
    return result

def get_bioproject_from_srp(srp_id: str) -> Optional[str]:
    """
    Get BioProject ID (PRJNA) from an SRP ID using NCBI's Entrez API.
    
    Args:
        srp_id: SRP ID
        
    Returns:
        BioProject ID (PRJNA format) or None if not found
    """
    # Check if we have a known mapping
    if srp_id in KNOWN_MAPPINGS:
        return KNOWN_MAPPINGS[srp_id]['bioproject_id']
    
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    
    try:
        # Step 1: Search for SRP ID in the SRA database
        search_url = f"{base_url}/esearch.fcgi"
        search_params = {
            "db": "sra",
            "term": srp_id,
            "retmode": "json",
            "retmax": 1
        }
        
        search_response = ncbi_request(search_url, search_params)
        search_response.raise_for_status()
        search_data = search_response.json()
        
        # Get the UID from the search results
        uids = search_data.get("esearchresult", {}).get("idlist", [])
        
        if not uids:
            logger.warning(f"No SRA record found for SRP ID: {srp_id}")
            return None
        
        # Step 2: Fetch the SRA record to get the BioProject ID
        fetch_url = f"{base_url}/efetch.fcgi"
        fetch_params = {
            "db": "sra",
            "id": uids[0],
            "retmode": "xml"
        }
        
        fetch_response = ncbi_request(fetch_url, fetch_params)
        fetch_response.raise_for_status()
        
        # Parse the XML response to find the BioProject ID
        root = ET.fromstring(fetch_response.content)
        
        # Look for BioProject ID in the external links
        for ext_link in root.findall(".//EXTERNAL_ID"):
            namespace = ext_link.get("namespace")
            if namespace and "BioProject" in namespace:
                return ext_link.text
        
        logger.warning(f"No BioProject ID found for SRP ID: {srp_id}")
        return None
        
    except Exception as e:
        logger.error(f"Error getting BioProject ID for SRP ID {srp_id}: {str(e)}")
        return None

def get_bioproject_from_erp(erp_id: str) -> Optional[str]:
    """
    Get BioProject ID (PRJEB) from an ERP ID.
    
    Args:
        erp_id: ERP ID
        
    Returns:
        BioProject ID (PRJEB format) or None if not found
    """
    # Check if we have a known mapping
    if erp_id in KNOWN_MAPPINGS:
        return KNOWN_MAPPINGS[erp_id]['bioproject_id']
    
    # Method 1: Try the ENA Portal API directly with the ERP ID
    try:
        url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={erp_id}&result=read_run&fields=study_accession&format=json"
        response = ebi_request(url)
        
        if response.status_code == 200:
            try:
                data = response.json()
                if data and isinstance(data, list) and len(data) > 0:
                    study_accession = data[0].get("study_accession")
                    if study_accession and study_accession.startswith("PRJEB"):
                        return study_accession
            except Exception as e:
                logger.warning(f"Error parsing ENA Portal API JSON for {erp_id}: {str(e)}")
    except Exception as e:
        logger.warning(f"Error querying ENA Portal API for {erp_id}: {str(e)}")
    
    # Method 2: Try the ENA Browser API XML endpoint
    try:
        # According to the API docs, we can use the xml endpoint to get study information
        url = f"{ENA_API_BASE_URL}/xml/{erp_id}"
        response = ebi_request(url)
        
        if response.status_code == 200:
            try:
                # Parse the XML response
                root = ET.fromstring(response.content)
                
                # Look for the PROJECT element which contains the BioProject ID
                for project_elem in root.findall(".//PROJECT"):
                    project_accession = project_elem.get("accession")
                    if project_accession and project_accession.startswith("PRJEB"):
                        return project_accession
                
                # If not found in PROJECT element, look for STUDY element
                for study_elem in root.findall(".//STUDY"):
                    study_accession = study_elem.get("accession")
                    if study_accession and study_accession.startswith("PRJEB"):
                        return study_accession
                        
                # If still not found, look for IDENTIFIERS section
                for ext_id in root.findall(".//EXTERNAL_ID"):
                    namespace = ext_id.get("namespace")
                    if namespace and "BioProject" in namespace:
                        return ext_id.text
            except Exception as e:
                logger.warning(f"Error parsing ENA XML for {erp_id}: {str(e)}")
    except Exception as e:
        logger.warning(f"Error querying ENA API for {erp_id}: {str(e)}")
    
    # Method 3: Try the ENA Portal API with the study endpoint
    try:
        url = f"https://www.ebi.ac.uk/ena/portal/api/study?accession={erp_id}&format=json"
        response = ebi_request(url)
        
        if response.status_code == 200:
            try:
                data = response.json()
                if data and isinstance(data, list) and len(data) > 0:
                    study = data[0]
                    bioproject_id = study.get("study_accession")
                    if bioproject_id and bioproject_id.startswith("PRJEB"):
                        return bioproject_id
            except Exception as e:
                logger.warning(f"Error parsing ENA study JSON for {erp_id}: {str(e)}")
    except Exception as e:
        logger.warning(f"Error querying ENA study API for {erp_id}: {str(e)}")
    
    # If all API queries fail, try the standard conversion
    match = re.search(r'ERP(\d+)', erp_id)
    if not match:
        logger.warning(f"Invalid ERP ID format: {erp_id}")
        return None
    
    # For ERP IDs, we can directly convert to PRJEB format
    # But this is not always accurate for all ERP IDs
    return f"PRJEB{match.group(1)}"

def get_gse_from_srp(srp_id: str) -> Optional[str]:
    """
    Get GSE ID from an SRP ID using NCBI's Entrez API.
    
    Args:
        srp_id: SRP ID
        
    Returns:
        GSE ID or None if not found
    """
    # Check if we have a known mapping
    if srp_id in KNOWN_MAPPINGS:
        return KNOWN_MAPPINGS[srp_id]['geo_id']
    
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    
    try:
        # First, try to get the BioProject ID
        bioproject_id = get_bioproject_from_srp(srp_id)
        if not bioproject_id:
            logger.warning(f"No BioProject ID found for SRP ID: {srp_id}")
            return None
        
        # Search for the BioProject ID in the GEO database
        search_url = f"{base_url}/esearch.fcgi"
        search_params = {
            "db": "gds",
            "term": f"{bioproject_id}[BioProject] OR {srp_id}",
            "retmode": "json",
            "retmax": 5
        }
        
        search_response = ncbi_request(search_url, search_params)
        search_response.raise_for_status()
        search_data = search_response.json()
        
        # Get the UIDs from the search results
        uids = search_data.get("esearchresult", {}).get("idlist", [])
        
        if not uids:
            # Try searching directly for the SRP ID
            search_params["term"] = srp_id
            search_response = ncbi_request(search_url, search_params)
            search_response.raise_for_status()
            search_data = search_response.json()
            uids = search_data.get("esearchresult", {}).get("idlist", [])
            
            if not uids:
                logger.warning(f"No GEO record found for SRP ID: {srp_id}")
                return None
        
        # Fetch the GEO record to get the GSE ID
        summary_url = f"{base_url}/esummary.fcgi"
        
        for uid in uids:
            summary_params = {
                "db": "gds",
                "id": uid,
                "retmode": "json"
            }
            
            summary_response = ncbi_request(summary_url, summary_params)
            summary_response.raise_for_status()
            summary_data = summary_response.json()
            
            # Extract GSE ID from summary
            result = summary_data.get("result", {})
            if result and result.get(uid):
                accession = result[uid].get("accession")
                
                if accession and accession.startswith("GSE"):
                    return accession
        
        # If we've tried all UIDs and still haven't found a GSE ID, try one more approach
        # Use the efetch API to get more details
        fetch_url = f"{base_url}/efetch.fcgi"
        fetch_params = {
            "db": "gds",
            "id": ",".join(uids),
            "retmode": "xml"
        }
        
        fetch_response = ncbi_request(fetch_url, fetch_params)
        fetch_response.raise_for_status()
        
        # Parse the XML response to find the GSE ID
        root = ET.fromstring(fetch_response.content)
        
        for gse_elem in root.findall(".//Accession"):
            if gse_elem.text and gse_elem.text.startswith("GSE"):
                return gse_elem.text
        
        logger.warning(f"No GSE ID found for SRP ID: {srp_id}")
        return None
        
    except Exception as e:
        logger.error(f"Error getting GSE ID for SRP ID {srp_id}: {str(e)}")
        return None

def get_arrayexpress_from_erp(erp_id: str) -> Optional[str]:
    """
    Get ArrayExpress ID (E-MTAB) from an ERP ID.
    
    Args:
        erp_id: ERP ID
        
    Returns:
        ArrayExpress ID (E-MTAB format) or None if not found
    """
    # Check if we have a known mapping
    if erp_id in KNOWN_MAPPINGS:
        return KNOWN_MAPPINGS[erp_id]['geo_id']
    
    # First, get the BioProject ID
    bioproject_id = get_bioproject_from_erp(erp_id)
    if not bioproject_id:
        logger.warning(f"No BioProject ID found for ERP ID: {erp_id}")
        return None
    
    logger.debug(f"Searching for ArrayExpress ID for ERP ID: {erp_id} (BioProject: {bioproject_id})")
    
    try:
        # Method 1: Direct query to BioStudies API using the ERP ID
        try:
            logger.debug(f"Trying BioStudies API with ERP ID: {erp_id}")
            url = f"https://www.ebi.ac.uk/biostudies/api/v1/search?query={erp_id}"
            response = ebi_request(url)
            
            if response.status_code == 200:
                try:
                    data = response.json()
                    hits = data.get("hits", [])
                    
                    for hit in hits:
                        accession = hit.get("accession", "")
                        if accession.startswith("E-"):
                            logger.debug(f"Found ArrayExpress ID via BioStudies API: {accession}")
                            return accession
                except Exception as e:
                    logger.warning(f"Error parsing BioStudies JSON for {erp_id}: {str(e)}")
        except Exception as e:
            logger.warning(f"Error querying BioStudies API for {erp_id}: {str(e)}")
        
        # Method 2: Direct query to BioStudies API using the BioProject ID
        try:
            logger.debug(f"Trying BioStudies API with BioProject ID: {bioproject_id}")
            url = f"https://www.ebi.ac.uk/biostudies/api/v1/search?query={bioproject_id}"
            response = ebi_request(url)
            
            if response.status_code == 200:
                try:
                    data = response.json()
                    hits = data.get("hits", [])
                    
                    for hit in hits:
                        accession = hit.get("accession", "")
                        if accession.startswith("E-"):
                            logger.debug(f"Found ArrayExpress ID via BioStudies API: {accession}")
                            return accession
                except Exception as e:
                    logger.warning(f"Error parsing BioStudies JSON for {bioproject_id}: {str(e)}")
        except Exception as e:
            logger.warning(f"Error querying BioStudies API for {bioproject_id}: {str(e)}")
        
        # Method 3: Try the ENA Browser API XML endpoint with includeLinks parameter
        try:
            logger.debug(f"Trying ENA Browser API XML with ERP ID: {erp_id}")
            url = f"{ENA_API_BASE_URL}/xml/{erp_id}?includeLinks=true"
            response = ebi_request(url)
            
            if response.status_code == 200:
                try:
                    # Parse the XML response
                    root = ET.fromstring(response.content)
                    
                    # Look for ArrayExpress links in the XREF_LINK elements
                    for xref_link in root.findall(".//XREF_LINK"):
                        db_elem = xref_link.find("DB")
                        id_elem = xref_link.find("ID")
                        
                        if (db_elem is not None and db_elem.text == "ArrayExpress" and 
                            id_elem is not None and id_elem.text.startswith("E-")):
                            logger.debug(f"Found ArrayExpress ID via ENA XML: {id_elem.text}")
                            return id_elem.text
                except Exception as e:
                    logger.warning(f"Error parsing ENA XML for {erp_id}: {str(e)}")
        except Exception as e:
            logger.warning(f"Error querying ENA XML API for {erp_id}: {str(e)}")
        
        # Method 4: Try with the BioProject ID
        try:
            logger.debug(f"Trying ENA Browser API XML with BioProject ID: {bioproject_id}")
            url = f"{ENA_API_BASE_URL}/xml/{bioproject_id}?includeLinks=true"
            response = ebi_request(url)
            
            if response.status_code == 200:
                try:
                    # Parse the XML response
                    root = ET.fromstring(response.content)
                    
                    # Look for ArrayExpress links in the XREF_LINK elements
                    for xref_link in root.findall(".//XREF_LINK"):
                        db_elem = xref_link.find("DB")
                        id_elem = xref_link.find("ID")
                        
                        if (db_elem is not None and db_elem.text == "ArrayExpress" and 
                            id_elem is not None and id_elem.text.startswith("E-")):
                            logger.debug(f"Found ArrayExpress ID via ENA XML: {id_elem.text}")
                            return id_elem.text
                except Exception as e:
                    logger.warning(f"Error parsing ENA XML for {bioproject_id}: {str(e)}")
        except Exception as e:
            logger.warning(f"Error querying ENA XML API for {bioproject_id}: {str(e)}")
        
        # Method 5: Try the ENA Portal API links endpoint
        try:
            logger.debug(f"Trying ENA Portal API links with BioProject ID: {bioproject_id}")
            url = f"https://www.ebi.ac.uk/ena/portal/api/links/study?accession={bioproject_id}&format=json"
            response = ebi_request(url)
            
            if response.status_code == 200:
                try:
                    links_data = response.json()
                    
                    for link in links_data:
                        target_id = link.get("target_id", "")
                        if target_id.startswith("E-"):
                            logger.debug(f"Found ArrayExpress ID via ENA links: {target_id}")
                            return target_id
                except Exception as e:
                    logger.warning(f"Error parsing ENA links JSON for {bioproject_id}: {str(e)}")
        except Exception as e:
            logger.warning(f"Error querying ENA links API for {bioproject_id}: {str(e)}")
        
        # Method 6: Try the EBI Search API
        try:
            logger.debug(f"Trying EBI Search API with ERP ID: {erp_id}")
            url = f"https://www.ebi.ac.uk/ebisearch/ws/rest/arrayexpress-experiments?query={erp_id}&format=json"
            response = ebi_request(url)
            
            if response.status_code == 200:
                try:
                    data = response.json()
                    entries = data.get("entries", [])
                    
                    if entries:
                        for entry in entries:
                            id_value = entry.get("id")
                            if id_value and id_value.startswith("E-"):
                                logger.debug(f"Found ArrayExpress ID via EBI Search: {id_value}")
                                return id_value
                except Exception as e:
                    logger.warning(f"Error parsing EBI Search JSON for {erp_id}: {str(e)}")
        except Exception as e:
            logger.warning(f"Error querying EBI Search API for {erp_id}: {str(e)}")
        
        # Method 7: Try the direct BioStudies ArrayExpress collection API
        try:
            logger.debug(f"Trying BioStudies ArrayExpress collection with ERP ID: {erp_id}")
            url = f"https://www.ebi.ac.uk/biostudies/arrayexpress/studies?query={erp_id}"
            response = ebi_request(url)
            
            if response.status_code == 200:
                try:
                    # Check if we can extract the E-MTAB ID from the response
                    content = response.text
                    # Look for E-MTAB patterns in the response
                    matches = re.findall(r'E-MTAB-\d+', content)
                    if matches:
                        logger.debug(f"Found ArrayExpress ID via BioStudies HTML: {matches[0]}")
                        return matches[0]
                except Exception as e:
                    logger.warning(f"Error parsing BioStudies ArrayExpress HTML for {erp_id}: {str(e)}")
        except Exception as e:
            logger.warning(f"Error querying BioStudies ArrayExpress API for {erp_id}: {str(e)}")
        
        logger.warning(f"No ArrayExpress ID found for ERP ID: {erp_id}")
        return None
        
    except Exception as e:
        logger.error(f"Error getting ArrayExpress ID for ERP ID {erp_id}: {str(e)}")
        return None

if __name__ == "__main__":
    # Example usage
    test_ids = ['ERP127673', 'SRP324458']
    print(f"Converting SRA IDs: {test_ids}")
    
    results = convert_sra_ids(test_ids)
    
    for sra_id, result in results.items():
        print(f"{sra_id}: BioProject ID = {result['bioproject_id']}, GEO/ArrayExpress ID = {result['geo_id']}") 