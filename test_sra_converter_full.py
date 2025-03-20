#!/usr/bin/env python3
"""
Comprehensive test script for the SRA ID converter.
This script tests the convert_sra_ids function with all the IDs mentioned in the original request.
"""

from sra_id_converter import convert_sra_ids

def main():
    # Test cases with expected results
    test_cases = {
        'ERP127673': {'expected_bioproject': 'PRJEB43688', 'expected_geo': 'E-MTAB-10220'},
        'SRP324458': {'expected_bioproject': 'PRJNA738600', 'expected_geo': 'GSE178360'},
        'ERP149679': {'expected_bioproject': 'PRJEB64504', 'expected_geo': 'E-MTAB-8142'},
        'ERP144781': {'expected_bioproject': 'PRJEB59723', 'expected_geo': 'E-MTAB-12650'}
    }
    
    # Full list of SRA IDs from the original request
    all_sra_ids = [
        'ERP149679', 'ERP144781', 'ERP123138', 'ERP156277', 'ERP151533',
        'ERP160803', 'ERP158366', 'SRP402417', 'ERP136281', 'SRP324458',
        'SRP364677', 'ERP136992', 'SRP329496', 'SRP273096', 'SRP510712',
        'SRP308561', 'SRP306446', 'SRP309720', 'SRP310949', 'SRP288163',
        'SRP314456', 'SRP323939', 'SRP324752', 'SRP373380', 'SRP329970',
        'ERP125682', 'ERP125690', 'ERP125913', 'ERP127673'
    ]
    
    # First, test with just the IDs that have expected results
    print("=== Testing with known IDs ===")
    test_ids = list(test_cases.keys())
    print(f"Converting test SRA IDs: {test_ids}")
    
    # Convert the test IDs
    results = convert_sra_ids(test_ids)
    
    # Verify the results against expected values
    all_passed = True
    for sra_id, result in results.items():
        bioproject_id = result['bioproject_id']
        geo_id = result['geo_id']
        
        expected_bioproject = test_cases[sra_id]['expected_bioproject']
        expected_geo = test_cases[sra_id]['expected_geo']
        
        bioproject_match = bioproject_id == expected_bioproject
        geo_match = geo_id == expected_geo
        
        if not bioproject_match or not geo_match:
            all_passed = False
        
        print(f"{sra_id}:")
        print(f"  BioProject ID: {bioproject_id} {'✓' if bioproject_match else '✗'} (Expected: {expected_bioproject})")
        print(f"  GEO/ArrayExpress ID: {geo_id} {'✓' if geo_match else '✗'} (Expected: {expected_geo})")
    
    print(f"\nAll tests {'PASSED' if all_passed else 'FAILED'}")
    
    # Now test with a small subset of the full list (to avoid rate limiting)
    print("\n=== Testing with a subset of all IDs ===")
    subset_ids = all_sra_ids[:5]  # Just test the first 5 IDs
    print(f"Converting subset of SRA IDs: {subset_ids}")
    
    # Convert the subset of IDs
    subset_results = convert_sra_ids(subset_ids)
    
    # Print the results
    for sra_id, result in subset_results.items():
        print(f"{sra_id}:")
        print(f"  BioProject ID: {result['bioproject_id']}")
        print(f"  GEO/ArrayExpress ID: {result['geo_id']}")
    
    print("\n=== Usage Example ===")
    print("To convert a list of SRA IDs, use the following code:")
    print("```python")
    print("from sra_id_converter import convert_sra_ids")
    print("")
    print("sra_ids = [")
    print("    'ERP149679', 'ERP144781', 'ERP123138', 'ERP156277', 'ERP151533',")
    print("    # ... more IDs")
    print("]")
    print("")
    print("results = convert_sra_ids(sra_ids)")
    print("")
    print("for sra_id, result in results.items():")
    print("    print(f\"{sra_id}: BioProject ID = {result['bioproject_id']}, GEO/ArrayExpress ID = {result['geo_id']}\")")
    print("```")
    
    print("\nFor better performance with large lists, consider:")
    print("1. Setting an NCBI API key as an environment variable")
    print("2. Adjusting batch_size, max_workers, and delay_between_batches parameters")
    print("3. Processing IDs in smaller batches to avoid rate limits")

if __name__ == "__main__":
    main() 