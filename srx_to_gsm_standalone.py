#!/usr/bin/env python
"""
Standalone function to convert SRX to GSM without requiring the pysradb package.
"""

import os
import sys
import time
import json
import requests
import pandas as pd
from collections import OrderedDict
from json.decoder import JSONDecodeError
from typing import List, Optional, Union


def get_retmax(n_records, retmax=500):
    """Get retmax for eutils."""
    return range(0, n_records, retmax)


def _retry_response(base_url, payload, key, max_retries=10):
    """Retry API request with exponential backoff."""
    for i in range(max_retries):
        time.sleep(2 ** i)  # Exponential backoff
        request = requests.get(base_url, params=OrderedDict(payload))
        try:
            response = request.json()
            if key in response and "error" not in response:
                return response
        except JSONDecodeError:
            continue
    sys.stderr.write(f"Failed to get response after {max_retries} retries\n")
    sys.exit(1)


def batch_srx_to_gsm(srx_list: List[str], batch_size: int = 200, api_key: Optional[str] = None) -> pd.DataFrame:
    """
    Convert a large list of SRX accessions to GSM accessions in batches.
    
    Parameters
    ----------
    srx_list : List[str]
        List of SRX accessions
    batch_size : int, optional
        Number of SRX accessions to process in each batch (default: 200)
    api_key : str, optional
        NCBI API key for faster requests
        
    Returns
    -------
    pandas.DataFrame
        DataFrame with SRX and corresponding GSM accessions
    """
    # Deduplicate SRX IDs to avoid redundant queries
    unique_srx = list(set(srx_list))
    
    # Process in batches
    results_dfs = []
    
    for i in range(0, len(unique_srx), batch_size):
        batch = unique_srx[i:i+batch_size]
        print(f"Processing batch {i//batch_size + 1}/{(len(unique_srx) + batch_size - 1)//batch_size} "
              f"({len(batch)} SRX IDs)")
        
        # Process this batch
        batch_result = srx_to_gsm(batch, api_key=api_key)
        
        if batch_result is not None and not batch_result.empty:
            results_dfs.append(batch_result)
        
        # Be nice to the API - add a delay between batches
        if i + batch_size < len(unique_srx):
            delay = 3 if api_key is None else 1
            print(f"Waiting {delay} seconds before next batch...")
            time.sleep(delay)
    
    # Combine all results
    if not results_dfs:
        print("No results found for any SRX IDs")
        return pd.DataFrame(columns=["experiment_accession", "experiment_alias"])
    
    combined_df = pd.concat(results_dfs, ignore_index=True)
    
    # Map back to original order if needed
    result_dict = dict(zip(combined_df["experiment_accession"], combined_df["experiment_alias"]))
    
    # Create a DataFrame with all original SRX IDs, including those without matches
    final_df = pd.DataFrame({
        "experiment_accession": srx_list,
        "experiment_alias": [result_dict.get(srx, None) for srx in srx_list]
    })
    
    return final_df


def srx_to_gsm(srx, api_key=None):
    """
    Convert SRX accession(s) to GSM accession(s).
    
    Parameters
    ----------
    srx : str or list
        SRX accession or list of SRX accessions
    api_key : str, optional
        NCBI API key for faster requests
        
    Returns
    -------
    pandas.DataFrame
        DataFrame with SRX and corresponding GSM accessions
    """
    if isinstance(srx, str):
        srx = [srx]
    
    # Set up base URLs and parameters
    base_url_esearch = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    base_url_esummary = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    
    esearch_params = [
        ("db", "gds"),
        ("usehistory", "n"),
        ("retmode", "json"),
    ]
    
    # Add API key if provided
    if api_key is not None:
        esearch_params.append(("api_key", str(api_key)))
        sleep_time = 1 / 10
    else:
        sleep_time = 1 / 3
    
    # Convert list of SRXs to search term
    term = " OR ".join(srx)
    
    # Initial search request
    payload = esearch_params.copy()
    payload += [("term", term)]
    request = requests.post(base_url_esearch, data=OrderedDict(payload))
    
    try:
        esearch_response = request.json()
    except JSONDecodeError:
        sys.stderr.write(
            f"Unable to parse esearch response json: {request.text}{os.linesep}. Will retry once.\n"
        )
        retry_after = request.headers.get("Retry-After", 1)
        time.sleep(int(retry_after))
        request = requests.post(base_url_esearch, data=OrderedDict(payload))
        try:
            esearch_response = request.json()
        except JSONDecodeError:
            sys.stderr.write(
                f"Unable to parse esearch response json: {request.text}{os.linesep}. Aborting.\n"
            )
            return None
    
    if "esummaryresult" in esearch_response:
        print("No result found")
        return None
    
    if "error" in esearch_response:
        # API rate limit exceeded
        esearch_response = _retry_response(
            base_url_esearch, payload, "esearchresult"
        )
    
    # Get number of records
    n_records = int(esearch_response["esearchresult"]["count"])
    
    # Create parameters for summary request
    def create_esummary_params(esearchresult):
        query_key = esearchresult["querykey"]
        webenv = esearchresult["webenv"]
        retstart = esearchresult["retstart"]
        retmax = 500
        
        return [
            ("query_key", query_key),
            ("WebEnv", webenv),
            ("retstart", retstart),
            ("retmax", retmax),
        ]
    
    # Get summary results
    results = {}
    for retstart in get_retmax(n_records):
        payload = esearch_params.copy()
        payload += create_esummary_params(esearch_response["esearchresult"])
        payload = OrderedDict(payload)
        payload["retstart"] = retstart
        request = requests.get(base_url_esummary, params=OrderedDict(payload))
        
        try:
            response = request.json()
        except JSONDecodeError:
            time.sleep(1)
            response = _retry_response(base_url_esummary, payload, "result")
        
        if "error" in response:
            # API rate limit exceeded
            response = _retry_response(base_url_esummary, payload, "result")
        
        if retstart == 0:
            results = response["result"]
        else:
            result = response["result"]
            for key, value in result.items():
                if key in list(results.keys()):
                    results[key] += value
                else:
                    results[key] = value
    
    # Process results
    try:
        uids = results["uids"]
    except KeyError:
        print(f"No results found for {srx} | Obtained result: {results}")
        return None
    
    gse_records = []
    for uid in uids:
        record = results[uid]
        del record["uid"]
        if record["extrelations"]:
            extrelations = record["extrelations"]
            for extrelation in extrelations:
                keys = list(extrelation.keys())
                values = list(extrelation.values())
                assert sorted(keys) == sorted(
                    ["relationtype", "targetobject", "targetftplink"]
                )
                assert len(values) == 3
                record[extrelation["relationtype"]] = extrelation["targetobject"]
            del record["extrelations"]
            gse_records.append(record)
    
    if not len(gse_records):
        print(f"No results found for {srx}")
        return None
    
    # Create DataFrame and filter for GSM entries
    gsm_df = pd.DataFrame(gse_records)
    gsm_df = gsm_df[gsm_df.entrytype == "GSM"].rename(
        columns={"SRA": "experiment_accession", "accession": "experiment_alias"}
    )
    gsm_df = gsm_df.loc[gsm_df["experiment_accession"].isin(srx)]
    
    # Return only the relevant columns
    return gsm_df[["experiment_accession", "experiment_alias"]].drop_duplicates()


def main():
    """Command line interface."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Convert SRX to GSM accessions")
    parser.add_argument("srx", nargs="+", help="SRX accession(s)")
    parser.add_argument("--api-key", help="NCBI API key")
    parser.add_argument("--batch-size", type=int, default=200, 
                        help="Number of SRX accessions to process in each batch (default: 200)")
    parser.add_argument("--output", help="Output file path (CSV)")
    args = parser.parse_args()
    
    # Use batch processing if there are many SRX IDs
    if len(args.srx) > args.batch_size:
        result = batch_srx_to_gsm(args.srx, batch_size=args.batch_size, api_key=args.api_key)
    else:
        result = srx_to_gsm(args.srx, api_key=args.api_key)
    
    if result is not None and not result.empty:
        if args.output:
            result.to_csv(args.output, index=False)
            print(f"Results saved to {args.output}")
        else:
            print(result.to_csv(index=False))


if __name__ == "__main__":
    main() 