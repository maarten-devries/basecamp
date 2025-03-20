#!/usr/bin/env python3

from sra_id_converter import convert_sra_ids

def main():
    # Test with original test cases
    test_ids = ['ERP149679', 'ERP144781', 'SRP324458']
    print(f"Testing with original IDs: {test_ids}")
    
    results = convert_sra_ids(test_ids)
    
    print("\nResults for original IDs:")
    for sra_id, result in results.items():
        print(f"{sra_id}: BioProject ID = {result['bioproject_id']}, GEO/ArrayExpress ID = {result['geo_id']}")
    
    # Test with additional IDs
    additional_ids = ['ERP127673', 'ERP136992', 'SRP323939', 'ERP125913']
    print(f"\nTesting with additional IDs: {additional_ids}")
    
    additional_results = convert_sra_ids(additional_ids)
    
    print("\nResults for additional IDs:")
    for sra_id, result in additional_results.items():
        print(f"{sra_id}: BioProject ID = {result['bioproject_id']}, GEO/ArrayExpress ID = {result['geo_id']}")

if __name__ == "__main__":
    main() 