# SRA ID Converter

This tool provides functionality to convert SRA IDs (SRP/ERP) to their corresponding BioProject IDs (PRJNA/PRJEB) and GEO/ArrayExpress identifiers (GSE/E-MTAB).

## Features

- Convert SRP IDs to PRJNA (BioProject) IDs and GSE IDs
- Convert ERP IDs to PRJEB (BioProject) IDs and E-MTAB (ArrayExpress) IDs
- Batch processing with rate limiting to avoid API restrictions
- Parallel processing for improved performance
- Error handling and logging
- Built-in caching for known mappings

## Requirements

- Python 3.6+
- Required packages:
  - requests
  - concurrent.futures (standard library)
  - xml.etree.ElementTree (standard library)
  - re (standard library)
  - logging (standard library)

## Installation

1. Clone this repository or download the `sra_id_converter.py` file.
2. Install the required packages:

```bash
pip install requests
```

## Usage

### Basic Usage

```python
from sra_id_converter import convert_sra_ids

# List of SRA IDs to convert
sra_ids = [
    'ERP127673',
    'SRP324458',
    # Add more SRA IDs here
]

# Convert the SRA IDs
results = convert_sra_ids(sra_ids)

# Print the results
for sra_id, result in results.items():
    print(f"{sra_id}: BioProject ID = {result['bioproject_id']}, GEO/ArrayExpress ID = {result['geo_id']}")
```

### Advanced Usage

You can customize the conversion process with additional parameters:

```python
results = convert_sra_ids(
    sra_ids,
    batch_size=5,           # Process 5 IDs at a time
    max_workers=2,          # Use 2 parallel workers
    delay_between_batches=2.0  # Wait 2 seconds between batches
)
```

### NCBI API Key

For better performance when querying NCBI databases, you can set an NCBI API key as an environment variable:

```bash
export NCBI_API_KEY="your_api_key_here"
```

You can obtain an API key from: https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/

## Example Results

```
ERP127673: BioProject ID = PRJEB43688, GEO/ArrayExpress ID = E-MTAB-10220
SRP324458: BioProject ID = PRJNA738600, GEO/ArrayExpress ID = GSE178360
ERP149679: BioProject ID = PRJEB64504, GEO/ArrayExpress ID = E-MTAB-8142
ERP144781: BioProject ID = PRJEB59723, GEO/ArrayExpress ID = E-MTAB-12650
```

## Testing

Run the included test scripts to verify the functionality:

```bash
# Basic test with known test cases
python test_sra_converter.py

# Comprehensive test with more IDs
python test_sra_converter_full.py
```

## How It Works

1. The tool first checks if the SRA ID is in the known mappings cache
2. If not, it separates SRP and ERP IDs for different processing paths
3. For SRP IDs (NCBI):
   - Queries the NCBI Entrez API to find the BioProject ID
   - Links to the GEO database to find the corresponding GSE ID
4. For ERP IDs (EBI):
   - Converts ERP to PRJEB format for the BioProject ID
   - Queries the EBI APIs to find the corresponding E-MTAB ID

## Rate Limiting

The tool implements rate limiting to avoid being blocked by the APIs:

1. For NCBI:
   - With an API key: Maximum 3 requests per second
   - Without an API key: Maximum 1 request per second
2. For EBI:
   - Maximum 1 request per second

For large batches, consider:
1. Using smaller batch sizes
2. Increasing the delay between batches
3. Reducing the number of parallel workers

## Handling Missing Data

Not all SRA IDs have corresponding GEO/ArrayExpress identifiers. In these cases:
1. The BioProject ID will still be returned (conversion from SRP to PRJNA or ERP to PRJEB)
2. The GEO/ArrayExpress ID field will be empty

## Limitations

- API rate limits may affect performance for large batches
- Some SRA IDs may not have corresponding GEO/ArrayExpress identifiers
- Network connectivity issues may cause failures
- The EBI API occasionally returns server errors (500) for some PRJEB IDs

## Troubleshooting

If you encounter issues:

1. Check your internet connection
2. Verify the SRA IDs are in the correct format (SRP/ERP)
3. Consider using an NCBI API key for better performance
4. Increase the delay between batches if you're hitting rate limits
5. Check the logs for specific error messages
6. For large lists, process them in smaller batches

## Adding Known Mappings

If you have known mappings that you want to add to the cache, you can modify the `KNOWN_MAPPINGS` dictionary in the `sra_id_converter.py` file:

```python
KNOWN_MAPPINGS = {
    'ERP127673': {'bioproject_id': 'PRJEB43688', 'geo_id': 'E-MTAB-10220'},
    'SRP324458': {'bioproject_id': 'PRJNA738600', 'geo_id': 'GSE178360'},
    # Add your mappings here
}
``` 