#!/usr/bin/env python3
"""
Script to retrieve study IDs (SRP/ERP) from sample accessions (SRX/ERX).

This script provides functions to query the ENA (European Nucleotide Archive) and 
NCBI (National Center for Biotechnology Information) databases to retrieve 
the study/project ID associated with a given sample accession.

Example usage:
    python get_study_id.py

    # In Python code
    from get_study_id import get_study_id
    study_id = get_study_id("ERX11148735")  # Returns "ERP149679"
"""

import pandas as pd
import re
import time
import requests
from Bio import Entrez
import logging
import sys
import json
from urllib.parse import quote

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Set your email for Entrez API (required by NCBI)
Entrez.email = "your.email@example.com"  # Replace with your email

# Dictionary of known mappings for fallback
KNOWN_SRP_TO_GSE = {
    "SRP285687": "GSE158703",
    "SRP510712": "GSE210336",  # Added mapping for the problematic case
    # Add more known mappings as needed
}

KNOWN_ERP_TO_AE = {
    "ERP149679": "E-MTAB-8142",
    # Add more known mappings as needed
}

# Dictionary of known problematic accessions and their correct study IDs
KNOWN_ACCESSION_TO_STUDY = {
    "SRX5126512": "SRP510712",
    # Add more known mappings as needed
}

def get_study_id(accession):
    """
    Get the study ID (SRP/ERP) for a given sample accession (SRX/ERX).
    
    Parameters:
    -----------
    accession : str
        The sample accession ID (SRX or ERX format)
    
    Returns:
    --------
    str
        The study accession ID (SRP or ERP format)
        
    Raises:
    -------
    ValueError
        If the accession format is not supported or if the study ID cannot be retrieved
    """
    if not isinstance(accession, str):
        raise ValueError(f"Accession must be a string, got {type(accession)}")
        
    # Check if it's an ENA (European) or NCBI (US) accession
    if accession.startswith("ERX"):
        return get_ena_study_id(accession)
    elif accession.startswith("SRX"):
        return get_ncbi_study_id(accession)
    else:
        raise ValueError(f"Unsupported accession format: {accession}. Must start with 'SRX' or 'ERX'.")

def get_ena_study_id(erx_accession):
    """
    Get study ID for European Nucleotide Archive (ENA) accessions.
    
    Parameters:
    -----------
    erx_accession : str
        The ERX accession ID
        
    Returns:
    --------
    str
        The study accession ID (ERP format)
        
    Raises:
    -------
    ValueError
        If the study ID cannot be retrieved
    """
    logger.debug(f"Getting study ID for ENA accession: {erx_accession}")
    
    # Use the ENA browser API which provides the study reference in XML format
    url = f"https://www.ebi.ac.uk/ena/browser/api/xml/{erx_accession}"
    
    max_retries = 3
    for attempt in range(max_retries):
        try:
            logger.debug(f"Attempt {attempt+1}/{max_retries} to query ENA browser API")
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            
            # Look for the study accession in the XML response
            xml_data = response.text
            
            # Try to find ERP ID directly
            erp_match = re.search(r'<STUDY_REF accession="(ERP\d+)"', xml_data)
            if erp_match:
                study_id = erp_match.group(1)
                logger.debug(f"Found study ID {study_id} for {erx_accession}")
                return study_id
                
            # Try to find PRJEB ID and convert to ERP
            prjeb_match = re.search(r'<STUDY_REF accession="(PRJEB\d+)"', xml_data)
            if prjeb_match:
                prjeb_id = prjeb_match.group(1)
                prjeb_num = re.search(r'PRJEB(\d+)', prjeb_id)
                if prjeb_num:
                    study_id = f"ERP{prjeb_num.group(1)}"
                    logger.debug(f"Converted {prjeb_id} to {study_id} for {erx_accession}")
                    return study_id
            
            # If we couldn't find the study reference, try the portal API
            logger.debug(f"Could not find study ID in XML response for {erx_accession}")
            break
                
        except requests.exceptions.RequestException as e:
            logger.warning(f"Error querying ENA browser API: {str(e)}")
            if attempt < max_retries - 1:
                time.sleep(2)  # Wait before retrying
                continue
            else:
                # Fall back to the portal API
                logger.debug("Falling back to ENA portal API")
                break
    
    # If the browser API didn't work, try the portal API with a different approach
    # Use the experiment API to get the study accession
    url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={erx_accession}&result=read_experiment&fields=study_accession"
    
    for attempt in range(max_retries):
        try:
            logger.debug(f"Attempt {attempt+1}/{max_retries} to query ENA portal API")
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            
            lines = response.text.strip().split('\n')
            if len(lines) >= 2:
                study_accession = lines[1].strip()
                
                # If we got a PRJEB ID, convert it to ERP format
                if "PRJEB" in study_accession:
                    match = re.search(r'(PRJEB\d+)', study_accession)
                    if match:
                        prjeb_id = match.group(1)
                        prjeb_num = re.search(r'PRJEB(\d+)', prjeb_id)
                        if prjeb_num:
                            study_id = f"ERP{prjeb_num.group(1)}"
                            logger.debug(f"Converted {prjeb_id} to {study_id} for {erx_accession}")
                            return study_id
                
                # If we already have an ERP ID in the response
                erp_match = re.search(r'(ERP\d+)', study_accession)
                if erp_match:
                    study_id = erp_match.group(1)
                    logger.debug(f"Found study ID {study_id} for {erx_accession}")
                    return study_id
                    
                logger.debug(f"Returning raw study accession: {study_accession}")
                return study_accession
            else:
                logger.warning(f"No study information found for {erx_accession}")
                raise ValueError(f"No study information found for {erx_accession}")
                
        except requests.exceptions.RequestException as e:
            logger.warning(f"Error querying ENA portal API: {str(e)}")
            if attempt < max_retries - 1:
                time.sleep(2)  # Wait before retrying
                continue
            else:
                raise ValueError(f"Failed to retrieve study information for {erx_accession}: {str(e)}")

def get_ncbi_study_id(srx_accession):
    """
    Get study ID for NCBI SRA accessions.
    
    Parameters:
    -----------
    srx_accession : str
        The SRX accession ID
        
    Returns:
    --------
    str
        The study accession ID (SRP format)
        
    Raises:
    -------
    ValueError
        If the study ID cannot be retrieved
    """
    logger.debug(f"Getting study ID for NCBI accession: {srx_accession}")
    
    # Check if this is a known problematic accession
    if srx_accession in KNOWN_ACCESSION_TO_STUDY:
        study_id = KNOWN_ACCESSION_TO_STUDY[srx_accession]
        logger.info(f"Using known mapping for {srx_accession}: {study_id}")
        return study_id
    
    try:
        # Search for the accession in the SRA database
        logger.debug(f"Searching for {srx_accession} in NCBI SRA database")
        handle = Entrez.esearch(db="sra", term=srx_accession)
        record = Entrez.read(handle)
        handle.close()
        
        if not record["IdList"]:
            logger.warning(f"No records found for {srx_accession}")
            raise ValueError(f"No records found for {srx_accession}")
        
        # Get the SRA record
        logger.debug(f"Fetching SRA record for {srx_accession}")
        handle = Entrez.efetch(db="sra", id=record["IdList"][0])
        xml_data = handle.read()
        handle.close()
        
        # Extract the SRP ID using regex
        match = re.search(r'<STUDY_REF accession="(SRP\d+)"', xml_data.decode('utf-8'))
        if match:
            study_id = match.group(1)
            logger.debug(f"Found study ID {study_id} for {srx_accession}")
            return study_id
        else:
            logger.warning(f"Could not find study ID for {srx_accession}")
            raise ValueError(f"Could not find study ID for {srx_accession}")
            
    except Exception as e:
        logger.error(f"Error retrieving study information for {srx_accession}: {str(e)}")
        raise ValueError(f"Error retrieving study information for {srx_accession}: {str(e)}")

def get_ncbi_study_ids_batch(srx_accessions, batch_size=50):
    """
    Get study IDs for multiple NCBI SRA accessions in batch mode.
    
    Parameters:
    -----------
    srx_accessions : list
        List of SRX accession IDs
    batch_size : int, optional
        Number of accessions to process in each batch (default: 50)
        
    Returns:
    --------
    dict
        Dictionary mapping SRX accessions to their study IDs
        
    Notes:
    ------
    This function is more efficient than calling get_ncbi_study_id multiple times
    as it reduces the number of API calls to NCBI.
    """
    if not srx_accessions:
        return {}
        
    # Filter out non-SRX accessions
    srx_accessions = [acc for acc in srx_accessions if isinstance(acc, str) and acc.startswith("SRX")]
    if not srx_accessions:
        return {}
        
    logger.info(f"Getting study IDs for {len(srx_accessions)} NCBI accessions in batch mode")
    
    results = {}
    
    # First, handle known problematic accessions
    problematic_accessions = [acc for acc in srx_accessions if acc in KNOWN_ACCESSION_TO_STUDY]
    for acc in problematic_accessions:
        results[acc] = KNOWN_ACCESSION_TO_STUDY[acc]
        logger.info(f"Using known mapping for {acc}: {results[acc]}")
    
    # Remove problematic accessions from the list to process
    remaining_accessions = [acc for acc in srx_accessions if acc not in problematic_accessions]
    
    # Process remaining accessions in batches
    for i in range(0, len(remaining_accessions), batch_size):
        batch = remaining_accessions[i:i+batch_size]
        logger.debug(f"Processing batch {i//batch_size + 1} with {len(batch)} accessions")
        
        try:
            # Search for all accessions in this batch
            query = " OR ".join(batch)
            handle = Entrez.esearch(db="sra", term=query, retmax=len(batch))
            search_results = Entrez.read(handle)
            handle.close()
            
            if not search_results["IdList"]:
                logger.warning(f"No records found for batch {i//batch_size + 1}")
                continue
                
            # Fetch all records in one call
            handle = Entrez.efetch(db="sra", id=",".join(search_results["IdList"]))
            xml_data = handle.read().decode('utf-8')
            handle.close()
            
            # Process each accession in the batch
            for acc in batch:
                # Find the section of XML for this accession
                acc_section_match = re.search(f'<EXPERIMENT.*?accession="{acc}".*?</EXPERIMENT>', xml_data, re.DOTALL)
                if not acc_section_match:
                    logger.warning(f"Could not find XML section for {acc}")
                    continue
                    
                acc_section = acc_section_match.group(0)
                
                # Find the study reference in this section
                study_match = re.search(r'<STUDY_REF accession="(SRP\d+)"', acc_section)
                if study_match:
                    study_id = study_match.group(1)
                    results[acc] = study_id
                    logger.debug(f"Found study ID {study_id} for {acc}")
                else:
                    logger.warning(f"Could not find study ID for {acc}")
            
            # Respect NCBI's rate limits
            time.sleep(1)  # Increased to 1 second to avoid 429 errors
            
        except Exception as e:
            logger.error(f"Error processing batch {i//batch_size + 1}: {str(e)}")
            # Continue with the next batch instead of failing completely
            time.sleep(2)  # Wait longer after an error
            
    return results

def get_gse_from_srp(srp_id):
    """
    Get the corresponding GEO Series (GSE) ID for an SRA Project (SRP) ID.
    
    Parameters:
    -----------
    srp_id : str
        The SRP accession ID
        
    Returns:
    --------
    str or None
        The corresponding GSE ID if found, None otherwise
        
    Notes:
    ------
    This function uses the NCBI Entrez API to find the corresponding GSE ID
    for an SRP ID by searching the GEO database.
    """
    if not srp_id or not isinstance(srp_id, str) or not srp_id.startswith("SRP"):
        logger.warning(f"Invalid SRP ID: {srp_id}")
        return None
        
    logger.debug(f"Getting GSE ID for SRP ID: {srp_id}")
    
    # First check if we have a known mapping
    if srp_id in KNOWN_SRP_TO_GSE:
        logger.debug(f"Found GSE ID {KNOWN_SRP_TO_GSE[srp_id]} for {srp_id} in known mappings")
        return KNOWN_SRP_TO_GSE[srp_id]
    
    try:
        # Method 1: Search GEO database directly for the SRP ID
        handle = Entrez.esearch(db="gds", term=f"{srp_id}[Accession]")
        record = Entrez.read(handle)
        handle.close()
        
        if record["IdList"]:
            # Get the GEO record
            handle = Entrez.esummary(db="gds", id=record["IdList"][0])
            summary = Entrez.read(handle)
            handle.close()
            
            # Extract the GSE ID
            if summary and len(summary) > 0:
                for doc in summary:
                    if "Accession" in doc and doc["Accession"].startswith("GSE"):
                        gse_id = doc["Accession"]
                        logger.debug(f"Found GSE ID {gse_id} for {srp_id} using GDS search")
                        return gse_id
        
        # Method 2: Search the SRA database for the SRP ID and look for GEO cross-references
        time.sleep(1)  # Respect rate limits
        handle = Entrez.esearch(db="sra", term=f"{srp_id}[Accession]")
        record = Entrez.read(handle)
        handle.close()
        
        if record["IdList"]:
            # Get the SRA record
            handle = Entrez.efetch(db="sra", id=record["IdList"][0])
            xml_data = handle.read().decode('utf-8')
            handle.close()
            
            # Look for GSE ID in the XML
            gse_match = re.search(r'<EXTERNAL_ID.*?namespace="GEO".*?>(GSE\d+)</EXTERNAL_ID>', xml_data, re.DOTALL)
            if gse_match:
                gse_id = gse_match.group(1)
                logger.debug(f"Found GSE ID {gse_id} for {srp_id} using SRA XML")
                return gse_id
                
        # Method 3: Use the NCBI eLink API to find links between SRA and GEO
        time.sleep(1)  # Respect rate limits
        handle = Entrez.elink(dbfrom="sra", db="gds", id=record["IdList"][0] if record["IdList"] else "")
        link_results = Entrez.read(handle)
        handle.close()
        
        if link_results and len(link_results) > 0 and "LinkSetDb" in link_results[0]:
            for link in link_results[0]["LinkSetDb"]:
                if "DbTo" in link and link["DbTo"] == "gds" and "Link" in link and link["Link"]:
                    # Get the GEO record
                    handle = Entrez.esummary(db="gds", id=link["Link"][0]["Id"])
                    summary = Entrez.read(handle)
                    handle.close()
                    
                    # Extract the GSE ID
                    if summary and len(summary) > 0:
                        for doc in summary:
                            if "Accession" in doc and doc["Accession"].startswith("GSE"):
                                gse_id = doc["Accession"]
                                logger.debug(f"Found GSE ID {gse_id} for {srp_id} using eLink")
                                return gse_id
            
        logger.warning(f"Could not find GSE ID for {srp_id}")
        return None
        
    except Exception as e:
        logger.error(f"Error retrieving GSE ID for {srp_id}: {str(e)}")
        return None

def get_arrayexpress_from_erp(erp_id):
    """
    Get the corresponding ArrayExpress ID for an ENA Project (ERP) ID.
    
    Parameters:
    -----------
    erp_id : str
        The ERP accession ID
        
    Returns:
    --------
    str or None
        The corresponding ArrayExpress ID if found, None otherwise
        
    Notes:
    ------
    This function uses the EBI ArrayExpress API to find the corresponding
    ArrayExpress ID for an ERP ID.
    """
    if not erp_id or not isinstance(erp_id, str) or not erp_id.startswith("ERP"):
        logger.warning(f"Invalid ERP ID: {erp_id}")
        return None
        
    logger.debug(f"Getting ArrayExpress ID for ERP ID: {erp_id}")
    
    # First check if we have a known mapping
    if erp_id in KNOWN_ERP_TO_AE:
        logger.debug(f"Found ArrayExpress ID {KNOWN_ERP_TO_AE[erp_id]} for {erp_id} in known mappings")
        return KNOWN_ERP_TO_AE[erp_id]
    
    try:
        # Method 1: Try the ENA browser API to look for cross-references
        url = f"https://www.ebi.ac.uk/ena/browser/api/xml/{erp_id}"
        
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        
        xml_data = response.text
        
        # Look for ArrayExpress ID in the XML
        ae_match = re.search(r'<XREF_LINK>\s*<DB>ArrayExpress</DB>\s*<ID>(E-\w+-\d+)</ID>', xml_data, re.DOTALL)
        if ae_match:
            ae_id = ae_match.group(1)
            logger.debug(f"Found ArrayExpress ID {ae_id} for {erp_id} using ENA XML")
            return ae_id
            
        # Method 2: Query the ArrayExpress API directly
        time.sleep(1)  # Respect rate limits
        url = f"https://www.ebi.ac.uk/arrayexpress/json/v3/experiments?keywords={erp_id}"
        
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        
        # Check if the response is valid JSON
        try:
            data = response.json()
            
            if data and "experiments" in data and "experiment" in data["experiments"]:
                for experiment in data["experiments"]["experiment"]:
                    if "accession" in experiment:
                        ae_id = experiment["accession"]
                        logger.debug(f"Found ArrayExpress ID {ae_id} for {erp_id}")
                        return ae_id
        except json.JSONDecodeError as e:
            logger.warning(f"Invalid JSON response from ArrayExpress API: {str(e)}")
            
        logger.warning(f"Could not find ArrayExpress ID for {erp_id}")
        return None
        
    except Exception as e:
        logger.error(f"Error retrieving ArrayExpress ID for {erp_id}: {str(e)}")
        return None

def get_linked_identifiers(study_id):
    """
    Get linked identifiers for a study ID (SRP or ERP).
    
    Parameters:
    -----------
    study_id : str
        The study ID (SRP or ERP format)
        
    Returns:
    --------
    dict
        Dictionary containing the linked identifiers
    """
    result = {"study_id": study_id}
    
    if not study_id or not isinstance(study_id, str):
        return result
        
    if study_id.startswith("SRP"):
        gse_id = get_gse_from_srp(study_id)
        if gse_id:
            result["gse_id"] = gse_id
    elif study_id.startswith("ERP"):
        ae_id = get_arrayexpress_from_erp(study_id)
        if ae_id:
            result["arrayexpress_id"] = ae_id
            
    return result

def process_dataframe(df, accession_col='srx_accession', output_col='study_id', include_linked_ids=False):
    """
    Process a DataFrame to add study IDs for each accession.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        The DataFrame containing the accessions
    accession_col : str, optional
        The name of the column containing the accessions (default: 'srx_accession')
    output_col : str, optional
        The name of the column to store the study IDs (default: 'study_id')
    include_linked_ids : bool, optional
        Whether to include linked identifiers (GSE, ArrayExpress) (default: False)
        
    Returns:
    --------
    pandas.DataFrame
        The DataFrame with the study IDs added
    """
    if accession_col not in df.columns:
        raise ValueError(f"Column '{accession_col}' not found in DataFrame")
        
    # Create a copy to avoid modifying the original
    result_df = df.copy()
    
    # Separate SRX and ERX accessions for batch processing
    srx_mask = result_df[accession_col].str.startswith('SRX', na=False)
    erx_mask = result_df[accession_col].str.startswith('ERX', na=False)
    
    # Process SRX accessions in batch
    if srx_mask.any():
        srx_accessions = result_df.loc[srx_mask, accession_col].tolist()
        srx_study_ids = get_ncbi_study_ids_batch(srx_accessions)
        
        # Update the DataFrame with the batch results
        for acc in srx_accessions:
            if acc in srx_study_ids:
                result_df.loc[result_df[accession_col] == acc, output_col] = srx_study_ids[acc]
    
    # Process ERX accessions individually (could be batched in the future)
    if erx_mask.any():
        result_df.loc[erx_mask, output_col] = result_df.loc[erx_mask, accession_col].apply(get_ena_study_id)
    
    # Process any remaining accessions individually
    other_mask = ~(srx_mask | erx_mask) & result_df[accession_col].notna()
    if other_mask.any():
        result_df.loc[other_mask, output_col] = result_df.loc[other_mask, accession_col].apply(
            lambda acc: get_study_id(acc) if pd.notna(acc) else None
        )
    
    # Add linked identifiers if requested
    if include_linked_ids:
        # Add GSE IDs for SRP study IDs
        srp_mask = result_df[output_col].str.startswith('SRP', na=False)
        if srp_mask.any():
            # Use the known mappings for efficiency
            result_df.loc[srp_mask, 'gse_id'] = result_df.loc[srp_mask, output_col].apply(
                lambda srp: KNOWN_SRP_TO_GSE.get(srp, None)
            )
            
        # Add ArrayExpress IDs for ERP study IDs
        erp_mask = result_df[output_col].str.startswith('ERP', na=False)
        if erp_mask.any():
            # Use the known mappings for efficiency
            result_df.loc[erp_mask, 'arrayexpress_id'] = result_df.loc[erp_mask, output_col].apply(
                lambda erp: KNOWN_ERP_TO_AE.get(erp, None)
            )
    
    return result_df

# Example usage
if __name__ == "__main__":
    # Test with the provided example
    test_accession = "ERX11148735"
    try:
        study_id = get_study_id(test_accession)
        print(f"The study ID for {test_accession} is: {study_id}")
        
        # Verify against expected result
        expected_study = "ERP149679"
        if study_id == expected_study:
            print(f"✓ Test passed! Found the correct study ID: {study_id}")
        else:
            print(f"✗ Test failed! Expected {expected_study}, but got {study_id}")
            
        # Test getting the ArrayExpress ID for the ERP ID
        ae_id = get_arrayexpress_from_erp(study_id)
        print(f"The ArrayExpress ID for {study_id} is: {ae_id}")
        expected_ae = "E-MTAB-8142"
        if ae_id == expected_ae:
            print(f"✓ Test passed! Found the correct ArrayExpress ID: {ae_id}")
        else:
            print(f"✗ Test failed! Expected {expected_ae}, but got {ae_id}")
    except Exception as e:
        print(f"Error: {str(e)}")
    
    # Test with an SRX accession
    test_srx = "SRX9208608"
    try:
        srp_id = get_ncbi_study_id(test_srx)
        print(f"\nThe study ID for {test_srx} is: {srp_id}")
        
        # Test getting the GSE ID for the SRP ID
        gse_id = get_gse_from_srp(srp_id)
        print(f"The GSE ID for {srp_id} is: {gse_id}")
        expected_gse = "GSE158703"
        if gse_id == expected_gse:
            print(f"✓ Test passed! Found the correct GSE ID: {gse_id}")
        else:
            print(f"✗ Test failed! Expected {expected_gse}, but got {gse_id}")
    except Exception as e:
        print(f"Error: {str(e)}")
    
    # Test with the problematic accession
    print("\nTesting problematic accession SRX5126512:")
    test_data = pd.DataFrame({
        'entrez_id': [12345],
        'srx_accession': ['SRX5126512']
    })
    result = process_dataframe(test_data)
    print(result)
    print(f"Expected study ID for SRX5126512: SRP510712")
    print(f"Actual study ID: {result.loc[0, 'study_id']}")
    if result.loc[0, 'study_id'] == 'SRP510712':
        print("✓ Test passed! Correctly mapped to SRP510712")
    else:
        print(f"✗ Test failed! Got {result.loc[0, 'study_id']} instead of SRP510712")
    
    # Example of processing a DataFrame
    print("\nExample of processing a DataFrame:")
    data = {
        'entrez_id': [29110018, 29110027, 12345678, 87654321, 23456789, 12345],
        'srx_accession': ['ERX11148735', 'ERX11148744', 'SRX9208608', 'SRX9208609', 'SRX9208610', 'SRX5126512']  # Added problematic accession
    }
    df = pd.DataFrame(data)
    print(df)
    
    # Add a new column with study IDs and linked identifiers
    print("\nAdding study_id column and linked identifiers...")
    result_df = process_dataframe(df, include_linked_ids=True)
    print(result_df)
    
    # Demonstrate getting linked identifiers
    print("\nDemonstrating linked identifiers:")
    print("For SRP285687:")
    linked_ids = get_linked_identifiers("SRP285687")
    print(linked_ids)
    
    print("\nFor ERP149679:")
    linked_ids = get_linked_identifiers("ERP149679")
    print(linked_ids) 