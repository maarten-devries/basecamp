from Bio import Entrez
import time
import re

def get_srp_for_srx_batch(srx_ids, email, batch_size=200, debug=False):
    """
    Get parent SRP IDs for a list of SRX IDs in batch mode.
    
    Args:
        srx_ids (list): List of SRX IDs
        email (str): Your email for Entrez API
        batch_size (int): Number of IDs to process in each batch
        debug (bool): Print debug information
        
    Returns:
        dict: Dictionary mapping SRX IDs to their parent SRP IDs
    """
    Entrez.email = email
    srx_to_srp = {}
    
    for i in range(0, len(srx_ids), batch_size):
        batch = srx_ids[i:i+batch_size]
        
        if debug:
            print(f"Processing batch: {batch}")
        
        # Use efetch to get the full record
        try:
            handle = Entrez.efetch(db="sra", id=",".join(batch), retmode="xml")
            xml_data = handle.read().decode('utf-8')
            handle.close()
            
            if debug:
                print(f"Received XML data of length: {len(xml_data)}")
            
            # Extract SRP IDs using regex
            for srx_id in batch:
                # Find the experiment section for this SRX
                exp_pattern = rf'<EXPERIMENT[^>]*accession="{srx_id}".*?</EXPERIMENT>'
                exp_match = re.search(exp_pattern, xml_data, re.DOTALL)
                
                if not exp_match:
                    if debug:
                        print(f"No experiment section found for {srx_id}")
                    continue
                
                exp_section = exp_match.group(0)
                
                # Find the study reference
                study_pattern = r'<STUDY_REF[^>]*accession="([^"]+)"'
                study_match = re.search(study_pattern, exp_section)
                
                if study_match:
                    srp_id = study_match.group(1)
                    srx_to_srp[srx_id] = srp_id
                    if debug:
                        print(f"Found mapping: {srx_id} -> {srp_id}")
                elif debug:
                    print(f"No study reference found for {srx_id}")
            
        except Exception as e:
            if debug:
                print(f"Error in efetch: {str(e)}")
        
        # Sleep to avoid overloading the API
        time.sleep(1)
    
    return srx_to_srp

# Example usage
if __name__ == "__main__":
    # Example with a single SRX ID
    srx_ids = ["SRX27443010"]
    result = get_srp_for_srx_batch(srx_ids, "your.email@example.com", debug=True)
    print(f"Result: {result}")
    
    # Example with multiple SRX IDs
    srx_ids_multiple = ["SRX27443010", "SRX27443011", "SRX27443012"]
    result_multiple = get_srp_for_srx_batch(srx_ids_multiple, "your.email@example.com")
    print(f"Result for multiple IDs: {result_multiple}") 