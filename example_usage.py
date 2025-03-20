import pandas as pd
from study_id_converter import add_geo_emtab_ids_to_dataframe
import os

# Example: Load your dataframe with study_ids
# In your case, you would load your actual dataframe instead
# result_df = pd.read_csv('your_data.csv')  # Uncomment and modify this line

# For demonstration, create a sample dataframe with some SRP and ERP IDs
data = {
    'study_id': [
        'SRP559437', 'SRP270870', 'ERP149679', 'SRP559437', 'SRP270870',
        'ERP149679', 'SRP559437', 'SRP270870', 'ERP149679', 'SRP559980'
    ]
}
result_df = pd.DataFrame(data)

print("Original dataframe:")
print(result_df)
print(f"Number of unique study IDs: {result_df['study_id'].nunique()}")

# Optional: Set NCBI API key as environment variable for better performance
# You can get an API key from: https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/
# os.environ['NCBI_API_KEY'] = 'your_api_key_here'

# Convert SRP/ERP IDs to GSE/E-MTAB IDs
# Using smaller batch size and longer delays to avoid rate limits
result_df = add_geo_emtab_ids_to_dataframe(
    df=result_df,
    study_id_col='study_id',
    batch_size=2,       # Process 2 unique IDs at a time to avoid rate limits
    max_workers=2,      # Use 2 parallel workers for API requests
    delay_between_batches=2.0  # Wait 2 seconds between batches
)

print("\nDataframe with added GSE/E-MTAB IDs:")
print(result_df)

# Show which SRP/ERP IDs were successfully converted
print("\nConversion results:")
for study_id, geo_id in zip(result_df['study_id'], result_df['geo_id']):
    if study_id != geo_id:
        print(f"{study_id} -> {geo_id}")
    else:
        print(f"{study_id} (not converted)")

# Count how many unique IDs were successfully converted
unique_pairs = result_df[['study_id', 'geo_id']].drop_duplicates()
converted = sum(unique_pairs['study_id'] != unique_pairs['geo_id'])
print(f"\nSuccessfully converted {converted} out of {unique_pairs.shape[0]} unique study IDs.")

# Save the results to a CSV file
result_df.to_csv('study_id_conversion_results.csv', index=False)
print("\nResults saved to study_id_conversion_results.csv") 