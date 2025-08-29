"""
UniProt Information Enhancer for Summary Data

This script enhances the existing summary.csv file by fetching additional protein information 
from UniProt including gene names, sequence lengths, intrinsically disordered region (IDR) 
lengths, and percentage of IDR content for each protein pair.

Usage:
    python extra_summary.py /path/to/base/directory

Requirements:
    - pandas library
    - requests library
    - Internet connection for UniProt API access
    - Existing summary.csv file in the base directory

Author: PLIA Project
"""

import pandas as pd
import requests
import time
import os
import argparse

def get_uniprot_info(uniprot_id):
    """
    Fetch protein information from UniProt REST API.
    
    Args:
        uniprot_id (str): UniProt accession ID (e.g., 'P12345')
    
    Returns:
        tuple: (gene_name, sequence_length, idr_total_length, percent_idr)
            - gene_name (str): Gene name or 'UNKNOWN' if not found
            - sequence_length (int): Total amino acid sequence length
            - idr_total_length (int): Total length of intrinsically disordered regions
            - percent_idr (float): Percentage of sequence that is disordered
    
    Note:
        Includes a 1-second delay between requests to be respectful to UniProt servers
    """
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    try:
        response = requests.get(url)
        response.raise_for_status()
        data = response.json()

        gene_name = "UNKNOWN"
        gene_section = data.get('genes', [])
        if gene_section and 'geneName' in gene_section[0]:
            gene_name = gene_section[0]['geneName']['value']

        sequence_length = data.get('sequence', {}).get('length', 0)

        idr_total_length = 0
        features = data.get('features', [])
        for feature in features:
            if feature.get('type') == 'Region' and feature.get('description') == 'Disordered':
                start = int(feature['location']['start']['value'])
                end = int(feature['location']['end']['value'])
                idr_total_length += end - start + 1

        percent_idr = (idr_total_length / sequence_length * 100) if sequence_length > 0 else 0

        return gene_name, sequence_length, idr_total_length, percent_idr

    except Exception as e:
        print(f"Error retrieving UniProt info for {uniprot_id}: {e}")
        return "UNKNOWN", None, None, None

def main(base_dir):
    """
    Main function that processes the summary.csv file and enhances it with UniProt data.
    
    This function:
    1. Reads the existing summary.csv file
    2. Extracts UniProt IDs from subfolder names
    3. Fetches protein information for each UniProt ID
    4. Adds new columns with gene names, sequence lengths, IDR information
    5. Saves the enhanced data back to summary.csv
    
    Args:
        base_dir (str): Path to directory containing summary.csv file
        
    Note:
        The function expects subfolder names to follow the pattern 'UNIPROT1_UNIPROT2_output'
        and will skip malformed entries while logging warnings.
    """
    input_path = os.path.join(base_dir, "summary.csv")

    df = pd.read_csv(input_path)

    gene_names_1, lengths_1, idrs_1, pidrs_1 = [], [], [], []
    gene_names_2, lengths_2, idrs_2, pidrs_2 = [], [], [], []
    gene_labels = []

    for subfolder in df['subfolder']:
        try:
            uniprot1, uniprot2 = subfolder.replace('_output', '').split('_')
        except ValueError:
            print(f"Skipping malformed subfolder: {subfolder}")
            gene_names_1.append("UNKNOWN")
            gene_names_2.append("UNKNOWN")
            lengths_1.append(None)
            lengths_2.append(None)
            idrs_1.append(None)
            idrs_2.append(None)
            pidrs_1.append(None)
            pidrs_2.append(None)
            gene_labels.append("UNKNOWN_UNKNOWN")
            continue

        gene1, len1, idr1, pidr1 = get_uniprot_info(uniprot1)
        gene2, len2, idr2, pidr2 = get_uniprot_info(uniprot2)

        gene_names_1.append(gene1)
        gene_names_2.append(gene2)
        lengths_1.append(len1)
        lengths_2.append(len2)
        idrs_1.append(idr1)
        idrs_2.append(idr2)
        pidrs_1.append(pidr1)
        pidrs_2.append(pidr2)
        gene_labels.append(f"{gene1}_{gene2}")

        time.sleep(1)  # Be respectful to UniProt servers

    # Add new columns to the dataframe
    df["gene_name_1"] = gene_names_1
    df["seq_len_1"] = lengths_1
    df["idr_len_1"] = idrs_1
    df["percent_idr_1"] = pidrs_1

    df["gene_name_2"] = gene_names_2
    df["seq_len_2"] = lengths_2
    df["idr_len_2"] = idrs_2
    df["percent_idr_2"] = pidrs_2

    df["gene_label"] = gene_labels

    output_path = os.path.join(base_dir, "summary.csv")
    df.to_csv(output_path, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process UniProt data from summary.csv")
    parser.add_argument("base_dir", help="Base directory containing summary.csv")
    args = parser.parse_args()
    main(args.base_dir)
