# GSE ID Finder

This repository contains scripts to find GSE IDs (GEO Series accession numbers) for single-cell RNA-seq data based on entrez_ids and SRX accessions.

## Requirements

- Python 3.6+
- Required packages:
  - pandas
  - requests
  - joblib (optional, for parallel processing)

Install the required packages with:

```bash
pip install pandas requests joblib
```

## Files

- `get_gse_id.py`: Contains functions to find GSE IDs using entrez_ids or SRX accessions
- `add_gse_ids.py`: Contains functions to add GSE IDs to a dataframe
- `process_csv.py`: Command-line script to process a CSV file and add GSE IDs

## Usage

### Basic Usage

To add GSE IDs to a CSV file containing entrez_ids and SRX accessions:

```bash
python process_csv.py input_file.csv
```

This will create a new file called `input_file_with_gse_ids.csv` with an additional `gse_id` column.

### Advanced Usage

```bash
python process_csv.py input_file.csv --output_file output_file.csv --entrez_id_col entrez_id --srx_col srx_accession --batch_size 20 --parallel --n_jobs 8
```

Options:
- `--output_file`: Specify the output file name (default: input_file_with_gse_ids.csv)
- `--entrez_id_col`: Specify the column name containing entrez_ids (default: entrez_id)
- `--srx_col`: Specify the column name containing SRX accessions (default: srx_accession)
- `--batch_size`: Number of rows to process before saving intermediate results (default: 10)
- `--parallel`: Use parallel processing
- `--n_jobs`: Number of parallel jobs (default: 4)

### Using the Python API

You can also use the functions directly in your Python code:

```python
import pandas as pd
from add_gse_ids import add_gse_ids_to_df

# Load your dataframe
df = pd.DataFrame({
    'entrez_id': [29110018, 29110027],
    'srx_accession': ['ERX11148735', 'ERX11148744']
})

# Add GSE IDs
result_df = add_gse_ids_to_df(df)

# Print the result
print(result_df)
```

## How It Works

The scripts use the NCBI Entrez API to query the SRA database and retrieve the GSE IDs associated with the entrez_ids or SRX accessions. The process works as follows:

1. For each row in the dataframe, the script first tries to find the GSE ID using the SRX accession (more reliable)
2. If that fails, it tries using the entrez_id
3. The GSE ID is added to the dataframe in a new column called `gse_id`

## Notes

- The NCBI Entrez API has rate limits, so the scripts include delays between requests
- For large datasets, the process can take a long time
- The parallel processing option can speed up the process, but be careful not to exceed the API rate limits
- The scripts save intermediate results to avoid losing progress if the process is interrupted 