#!/usr/bin/env python3
"""
Test script for the SRA ID converter.
This script demonstrates how to use the convert_sra_ids function
to retrieve BioProject IDs and GEO/ArrayExpress identifiers from SRA IDs.
"""

from sra_id_converter import convert_sra_ids

def main():
    # Test cases provided by the user
    test_cases = {
        'ERP127673': {'expected_bioproject': 'PRJEB43688', 'expected_geo': 'E-MTAB-10220'},
        'SRP324458': {'expected_bioproject': 'PRJNA738600', 'expected_geo': 'GSE178360'}
    }
    
    # Example list of SRA IDs
    sra_ids = [
        'ERP149679', 'ERP144781', 'ERP123138', 'ERP156277', 'ERP151533',
        'ERP160803', 'ERP158366', 'SRP402417', 'ERP136281', 'SRP324458',
        'SRP364677', 'ERP136992', 'SRP329496', 'SRP273096', 'SRP510712',
        'SRP308561', 'SRP306446', 'SRP309720', 'SRP310949', 'SRP288163',
        'SRP314456', 'SRP323939', 'SRP324752', 'SRP373380', 'SRP329970',
        'ERP125682', 'ERP125690', 'ERP125913', 'ERP127673'
    ]
    
    # For testing, we'll just use the specific test cases
    test_ids = list(test_cases.keys())
    print(f"Converting test SRA IDs: {test_ids}")
    
    # Convert the test IDs
    results = convert_sra_ids(test_ids, batch_size=2, max_workers=2, delay_between_batches=1.0)
    
    # Verify the results against expected values
    for sra_id, result in results.items():
        bioproject_id = result['bioproject_id']
        geo_id = result['geo_id']
        
        expected_bioproject = test_cases[sra_id]['expected_bioproject']
        expected_geo = test_cases[sra_id]['expected_geo']
        
        bioproject_match = bioproject_id == expected_bioproject
        geo_match = geo_id == expected_geo
        
        print(f"{sra_id}:")
        print(f"  BioProject ID: {bioproject_id} {'✓' if bioproject_match else '✗'} (Expected: {expected_bioproject})")
        print(f"  GEO/ArrayExpress ID: {geo_id} {'✓' if geo_match else '✗'} (Expected: {expected_geo})")
    
    print("\n--- Full Example ---")
    print("To convert a larger list of SRA IDs, you would use:")
    print("results = convert_sra_ids(sra_ids)")
    print("This would return a dictionary with the following structure:")
    print("{")
    print("  'ERP127673': {'bioproject_id': 'PRJEB43688', 'geo_id': 'E-MTAB-10220'},")
    print("  'SRP324458': {'bioproject_id': 'PRJNA738600', 'geo_id': 'GSE178360'},")
    print("  # ... more results")
    print("}")

if __name__ == "__main__":
    main() 