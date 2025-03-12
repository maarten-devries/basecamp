import requests
import time
import re
import xml.etree.ElementTree as ET
from typing import Optional, Dict, Any, Union

def get_gse_id_from_entrez(entrez_id: Union[str, int]) -> Optional[str]:
    """
    Find the GSE ID for a given entrez_id by querying the NCBI Entrez API.
    
    Args:
        entrez_id: The entrez_id to search for
        
    Returns:
        The GSE ID if found, None otherwise
    """
    # Convert entrez_id to string if it's an integer
    entrez_id = str(entrez_id)
    
    # Base URL for NCBI Entrez API
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    
    # First, search for the entrez_id in the SRA database
    search_url = f"{base_url}esearch.fcgi"
    search_params = {
        "db": "sra",
        "term": entrez_id,
        "retmax": 1,
        "usehistory": "y"
    }
    
    try:
        search_response = requests.get(search_url, params=search_params)
        search_response.raise_for_status()
        
        # Parse the XML response
        search_root = ET.fromstring(search_response.content)
        
        # Get the ID from the search result
        id_elements = search_root.findall(".//Id")
        if not id_elements:
            print(f"No SRA record found for entrez_id: {entrez_id}")
            return None
        
        sra_id = id_elements[0].text
        
        # Use the SRA ID to get the full record
        fetch_url = f"{base_url}efetch.fcgi"
        fetch_params = {
            "db": "sra",
            "id": sra_id,
            "retmode": "xml"
        }
        
        # Add a small delay to avoid hitting NCBI rate limits
        time.sleep(0.5)
        
        fetch_response = requests.get(fetch_url, params=fetch_params)
        fetch_response.raise_for_status()
        
        # Parse the XML response
        fetch_content = fetch_response.content.decode('utf-8')
        
        # Look for GSE ID in the XML content
        # GSE IDs typically follow the pattern "GSE" followed by numbers
        gse_match = re.search(r'(GSE\d+)', fetch_content)
        
        if gse_match:
            return gse_match.group(1)
        else:
            # Try an alternative approach - look for external_id with namespace="GEO"
            fetch_root = ET.fromstring(fetch_response.content)
            external_ids = fetch_root.findall(".//EXTERNAL_ID")
            
            for ext_id in external_ids:
                namespace = ext_id.get("namespace")
                if namespace and "GEO" in namespace:
                    return ext_id.text
            
            print(f"No GSE ID found for entrez_id: {entrez_id}")
            return None
            
    except requests.exceptions.RequestException as e:
        print(f"Error querying NCBI API: {e}")
        return None
    except ET.ParseError as e:
        print(f"Error parsing XML response: {e}")
        return None

def get_gse_id_from_srx(srx_accession: str) -> Optional[str]:
    """
    Find the GSE ID for a given SRX accession by querying the NCBI Entrez API.
    
    Args:
        srx_accession: The SRX accession to search for
        
    Returns:
        The GSE ID if found, None otherwise
    """
    # Base URL for NCBI Entrez API
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    
    # First, search for the SRX accession in the SRA database
    search_url = f"{base_url}esearch.fcgi"
    search_params = {
        "db": "sra",
        "term": srx_accession,
        "retmax": 1,
        "usehistory": "y"
    }
    
    try:
        search_response = requests.get(search_url, params=search_params)
        search_response.raise_for_status()
        
        # Parse the XML response
        search_root = ET.fromstring(search_response.content)
        
        # Get the ID from the search result
        id_elements = search_root.findall(".//Id")
        if not id_elements:
            print(f"No SRA record found for SRX accession: {srx_accession}")
            return None
        
        sra_id = id_elements[0].text
        
        # Use the SRA ID to get the full record
        fetch_url = f"{base_url}efetch.fcgi"
        fetch_params = {
            "db": "sra",
            "id": sra_id,
            "retmode": "xml"
        }
        
        # Add a small delay to avoid hitting NCBI rate limits
        time.sleep(0.5)
        
        fetch_response = requests.get(fetch_url, params=fetch_params)
        fetch_response.raise_for_status()
        
        # Parse the XML response
        fetch_content = fetch_response.content.decode('utf-8')
        
        # Look for GSE ID in the XML content
        # GSE IDs typically follow the pattern "GSE" followed by numbers
        gse_match = re.search(r'(GSE\d+)', fetch_content)
        
        if gse_match:
            return gse_match.group(1)
        else:
            # Try an alternative approach - look for external_id with namespace="GEO"
            fetch_root = ET.fromstring(fetch_response.content)
            external_ids = fetch_root.findall(".//EXTERNAL_ID")
            
            for ext_id in external_ids:
                namespace = ext_id.get("namespace")
                if namespace and "GEO" in namespace:
                    return ext_id.text
            
            # Try another approach - look for study_ref with namespace="GEO"
            study_refs = fetch_root.findall(".//STUDY_REF")
            for study_ref in study_refs:
                identifiers = study_ref.findall(".//IDENTIFIERS/EXTERNAL_ID")
                for identifier in identifiers:
                    namespace = identifier.get("namespace")
                    if namespace and "GEO" in namespace:
                        return identifier.text
            
            print(f"No GSE ID found for SRX accession: {srx_accession}")
            return None
            
    except requests.exceptions.RequestException as e:
        print(f"Error querying NCBI API: {e}")
        return None
    except ET.ParseError as e:
        print(f"Error parsing XML response: {e}")
        return None

def get_gse_id_direct(srx_accession: str) -> Optional[str]:
    """
    Find the GSE ID for a given SRX accession by directly searching the GEO database.
    
    Args:
        srx_accession: The SRX accession to search for
        
    Returns:
        The GSE ID if found, None otherwise
    """
    # Base URL for NCBI Entrez API
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    
    # Search for the SRX accession in the GEO database
    search_url = f"{base_url}esearch.fcgi"
    search_params = {
        "db": "gds",  # GEO DataSets database
        "term": srx_accession,
        "retmax": 10
    }
    
    try:
        search_response = requests.get(search_url, params=search_params)
        search_response.raise_for_status()
        
        # Parse the XML response
        search_root = ET.fromstring(search_response.content)
        
        # Get the IDs from the search result
        id_elements = search_root.findall(".//Id")
        if not id_elements:
            print(f"No GEO record found for SRX accession: {srx_accession}")
            return None
        
        # Get all IDs
        geo_ids = [id_elem.text for id_elem in id_elements]
        
        # Fetch summary for each ID
        summary_url = f"{base_url}esummary.fcgi"
        summary_params = {
            "db": "gds",
            "id": ",".join(geo_ids),
            "retmode": "xml"
        }
        
        # Add a small delay to avoid hitting NCBI rate limits
        time.sleep(0.5)
        
        summary_response = requests.get(summary_url, params=summary_params)
        summary_response.raise_for_status()
        
        # Parse the XML response
        summary_root = ET.fromstring(summary_response.content)
        
        # Look for GSE IDs in the summaries
        for doc_sum in summary_root.findall(".//DocSum"):
            # Get the accession
            accession = None
            for item in doc_sum.findall(".//Item"):
                if item.get("Name") == "Accession":
                    accession = item.text
                    if accession.startswith("GSE"):
                        return accession
        
        print(f"No GSE ID found for SRX accession: {srx_accession}")
        return None
            
    except requests.exceptions.RequestException as e:
        print(f"Error querying NCBI API: {e}")
        return None
    except ET.ParseError as e:
        print(f"Error parsing XML response: {e}")
        return None

def get_gse_id_from_bioproject(srx_accession: str) -> Optional[str]:
    """
    Find the GSE ID for a given SRX accession by first finding the BioProject ID,
    then searching for the GSE ID associated with that BioProject.
    
    Args:
        srx_accession: The SRX accession to search for
        
    Returns:
        The GSE ID if found, None otherwise
    """
    # Base URL for NCBI Entrez API
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    
    # First, search for the SRX accession in the SRA database
    search_url = f"{base_url}esearch.fcgi"
    search_params = {
        "db": "sra",
        "term": srx_accession,
        "retmax": 1
    }
    
    try:
        search_response = requests.get(search_url, params=search_params)
        search_response.raise_for_status()
        
        # Parse the XML response
        search_root = ET.fromstring(search_response.content)
        
        # Get the ID from the search result
        id_elements = search_root.findall(".//Id")
        if not id_elements:
            print(f"No SRA record found for SRX accession: {srx_accession}")
            return None
        
        sra_id = id_elements[0].text
        
        # Use the SRA ID to get the full record
        fetch_url = f"{base_url}efetch.fcgi"
        fetch_params = {
            "db": "sra",
            "id": sra_id,
            "retmode": "xml"
        }
        
        # Add a small delay to avoid hitting NCBI rate limits
        time.sleep(0.5)
        
        fetch_response = requests.get(fetch_url, params=fetch_params)
        fetch_response.raise_for_status()
        
        # Parse the XML response to find the BioProject ID
        fetch_root = ET.fromstring(fetch_response.content)
        
        # Look for BioProject ID
        bioproject_id = None
        external_ids = fetch_root.findall(".//EXTERNAL_ID")
        
        for ext_id in external_ids:
            namespace = ext_id.get("namespace")
            if namespace and "BioProject" in namespace:
                bioproject_id = ext_id.text
                break
        
        if not bioproject_id:
            print(f"No BioProject ID found for SRX accession: {srx_accession}")
            return None
        
        # Now search for the BioProject ID in the GEO database
        bioproject_search_params = {
            "db": "gds",  # GEO DataSets database
            "term": bioproject_id,
            "retmax": 10
        }
        
        # Add a small delay to avoid hitting NCBI rate limits
        time.sleep(0.5)
        
        bioproject_search_response = requests.get(search_url, params=bioproject_search_params)
        bioproject_search_response.raise_for_status()
        
        # Parse the XML response
        bioproject_search_root = ET.fromstring(bioproject_search_response.content)
        
        # Get the IDs from the search result
        bioproject_id_elements = bioproject_search_root.findall(".//Id")
        if not bioproject_id_elements:
            print(f"No GEO record found for BioProject ID: {bioproject_id}")
            return None
        
        # Get all IDs
        geo_ids = [id_elem.text for id_elem in bioproject_id_elements]
        
        # Fetch summary for each ID
        summary_url = f"{base_url}esummary.fcgi"
        summary_params = {
            "db": "gds",
            "id": ",".join(geo_ids),
            "retmode": "xml"
        }
        
        # Add a small delay to avoid hitting NCBI rate limits
        time.sleep(0.5)
        
        summary_response = requests.get(summary_url, params=summary_params)
        summary_response.raise_for_status()
        
        # Parse the XML response
        summary_root = ET.fromstring(summary_response.content)
        
        # Look for GSE IDs in the summaries
        for doc_sum in summary_root.findall(".//DocSum"):
            # Get the accession
            accession = None
            for item in doc_sum.findall(".//Item"):
                if item.get("Name") == "Accession":
                    accession = item.text
                    if accession.startswith("GSE"):
                        return accession
        
        print(f"No GSE ID found for BioProject ID: {bioproject_id}")
        return None
            
    except requests.exceptions.RequestException as e:
        print(f"Error querying NCBI API: {e}")
        return None
    except ET.ParseError as e:
        print(f"Error parsing XML response: {e}")
        return None

def get_gse_id_from_ena(srx_accession: str) -> Optional[str]:
    """
    Find the GSE ID for a given SRX/ERX accession by querying the ENA API.
    
    Args:
        srx_accession: The SRX/ERX accession to search for
        
    Returns:
        The GSE ID if found, None otherwise
    """
    print(f"Searching for GSE ID for {srx_accession} using ENA API...")
    
    # ENA API URL
    ena_url = "https://www.ebi.ac.uk/ena/portal/api/filereport"
    
    # Parameters for the API request
    params = {
        "accession": srx_accession,
        "result": "read_experiment",
        "fields": "study_accession,secondary_study_accession,experiment_accession,experiment_title,study_title"
    }
    
    try:
        # Make the API request
        print(f"Making request to {ena_url} with params: {params}")
        response = requests.get(ena_url, params=params)
        response.raise_for_status()
        
        # Print the response for debugging
        print(f"Response status code: {response.status_code}")
        print(f"Response content: {response.text}")
        
        # Parse the response
        lines = response.text.strip().split('\n')
        if len(lines) < 2:
            print(f"No ENA record found for accession: {srx_accession}")
            return None
        
        # Get the header and data
        header = lines[0].split('\t')
        data = lines[1].split('\t')
        
        print(f"Header: {header}")
        print(f"Data: {data}")
        
        # Create a dictionary from the header and data
        record = dict(zip(header, data))
        
        print(f"Record: {record}")
        
        # Look for GSE ID in the study_title or secondary_study_accession
        study_title = record.get('study_title', '')
        secondary_study_accession = record.get('secondary_study_accession', '')
        
        print(f"Study title: {study_title}")
        print(f"Secondary study accession: {secondary_study_accession}")
        
        # Check if the secondary_study_accession is a GSE ID
        if secondary_study_accession.startswith('GSE'):
            print(f"Found GSE ID in secondary_study_accession: {secondary_study_accession}")
            return secondary_study_accession
        
        # Look for GSE ID in the study_title
        gse_match = re.search(r'(GSE\d+)', study_title)
        if gse_match:
            gse_id = gse_match.group(1)
            print(f"Found GSE ID in study_title: {gse_id}")
            return gse_id
        
        print(f"No GSE ID found for accession: {srx_accession}")
        return None
        
    except requests.exceptions.RequestException as e:
        print(f"Error querying ENA API: {e}")
        return None

def get_gse_id(entrez_id=None, srx_accession=None) -> Optional[str]:
    """
    Find the GSE ID using multiple methods.
    
    Args:
        entrez_id: The entrez_id to search for
        srx_accession: The SRX accession to search for
        
    Returns:
        The GSE ID if found, None otherwise
    """
    # Try all available methods
    if srx_accession:
        # Try ENA API first for ERX accessions
        if srx_accession.startswith('ERX'):
            gse_id = get_gse_id_from_ena(srx_accession)
            if gse_id:
                return gse_id
        
        # Try direct GEO search
        gse_id = get_gse_id_direct(srx_accession)
        if gse_id:
            return gse_id
        
        # Try BioProject method
        gse_id = get_gse_id_from_bioproject(srx_accession)
        if gse_id:
            return gse_id
        
        # Try SRX method
        gse_id = get_gse_id_from_srx(srx_accession)
        if gse_id:
            return gse_id
    
    if entrez_id:
        # Try entrez_id method
        gse_id = get_gse_id_from_entrez(entrez_id)
        if gse_id:
            return gse_id
    
    return None

# Test function with an example
if __name__ == "__main__":
    # Test with entrez_id
    test_entrez_id = 29110018  # Example from your dataframe
    test_srx = "ERX11148735"  # Example from your dataframe
    
    print(f"Testing with entrez_id: {test_entrez_id} and SRX accession: {test_srx}")
    
    # Try all methods
    gse_id = get_gse_id(entrez_id=test_entrez_id, srx_accession=test_srx)
    print(f"GSE ID: {gse_id}") 