# Study ID Converter

This tool converts SRP/ERP IDs to their corresponding GSE/E-MTAB IDs by querying public databases.

## Features

- Converts SRP IDs to GSE IDs using NCBI's Entrez API
- Converts ERP IDs to E-MTAB IDs using EBI's ArrayExpress API
- Efficiently processes IDs in batches to reduce API calls
- Uses parallel processing for faster conversion
- Handles large datasets with thousands of IDs

## Requirements

```
pandas
requests
```

## Usage

### Basic Usage

```python
import pandas as pd
from study_id_converter import add_geo_emtab_ids_to_dataframe

# Load your dataframe with study_ids
result_df = pd.DataFrame({'study_id': ['SRP559437', 'SRP270870', 'ERP149679']})

# Convert SRP/ERP IDs to GSE/E-MTAB IDs
result_df = add_geo_emtab_ids_to_dataframe(
    df=result_df,
    study_id_col='study_id',
    batch_size=50,  # Process 50 unique IDs at a time
    max_workers=5   # Use 5 parallel workers for API requests
)

# The result_df now has a new 'geo_id' column with the converted IDs
print(result_df)
```

### Advanced Usage

You can also use the lower-level functions directly:

```python
from study_id_converter import convert_study_ids

# List of SRP/ERP IDs
study_ids = ['SRP559437', 'SRP270870', 'ERP149679']

# Convert IDs
id_mapping = convert_study_ids(
    study_ids=study_ids,
    batch_size=50,
    max_workers=5
)

# id_mapping is a dictionary mapping original IDs to converted IDs
print(id_mapping)
```

## How It Works

1. The tool separates SRP and ERP IDs for processing with different APIs
2. For SRP IDs:
   - Queries NCBI's Entrez API to find corresponding GEO datasets
   - Extracts GSE IDs from the search results
3. For ERP IDs:
   - Queries EBI's ArrayExpress API to find corresponding experiments
   - Extracts E-MTAB IDs from the search results
4. Results are combined and returned as a dictionary or added to the dataframe

## Performance Considerations

- The tool processes IDs in batches to reduce the number of API calls
- It uses parallel processing to speed up the conversion
- For very large datasets (thousands of IDs), the conversion may take several minutes
- Be respectful of the API rate limits by using reasonable batch sizes and adding delays between batches

## Example

See `example_usage.py` for a complete example of how to use this tool with a dataframe.

## Limitations

- Not all SRP/ERP IDs have corresponding GSE/E-MTAB IDs
- The conversion relies on public APIs which may change or have rate limits
- Some conversions may fail due to network issues or API limitations 