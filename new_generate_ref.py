"""
Domain Reference File Generator

This script generates domain reference files for protein pairs by fetching sequences from UniProt
and extracting domain regions based on InterPro and UniProt binding site annotations. It creates
reference CSV files for each protein pair directory containing known binding domains.

Usage:
    python new_generate_ref.py --base_dir /path/to/base/directory

Requirements:
    - pandas library
    - requests library
    - Internet connection for UniProt API access
    - CSV files with domain annotations:
        - human_interpro_ppi_domains_consolidated_min3.csv
        - human_uniprot_ppi_sites_min3.csv

Author: PLIA Project
"""

import os
import pandas as pd
import requests
import argparse

UNIPROT_API = "https://rest.uniprot.org/uniprotkb/{}.fasta"

def fetch_uniprot_sequence(uniprot_id):
    """
    Fetch protein sequence from UniProt REST API.
    
    Args:
        uniprot_id (str): UniProt accession ID (e.g., 'P12345')
    
    Returns:
        str or None: Protein sequence string, or None if fetch fails
        
    Note:
        Removes the FASTA header line and returns only the sequence
    """
    try:
        response = requests.get(UNIPROT_API.format(uniprot_id))
        response.raise_for_status()
        lines = response.text.strip().split("\n")
        return "".join(lines[1:])  # Remove header line
    except Exception as e:
        print(f"[ERROR] Fetching {uniprot_id} failed: {e}")
        return None

def extract_domains(uniprot_id, domain_df):
    """
    Extract domain sequences for a given UniProt ID.
    
    Args:
        uniprot_id (str): UniProt accession ID
        domain_df (pandas.DataFrame): DataFrame containing domain annotations with columns:
                                     UniProt_AC, Start, End
    
    Returns:
        list or None: List of domain sequences, or None if no domains found
                     or sequence fetch fails
        
    Note:
        Uses 1-indexed coordinates from UniProt annotations
    """
    domains = domain_df[domain_df["UniProt_AC"] == uniprot_id]
    if domains.empty:
        return None

    sequence = fetch_uniprot_sequence(uniprot_id)
    if not sequence:
        return None

    seqs = []
    for _, row in domains.iterrows():
        start, end = row["Start"], row["End"]
        if 1 <= start <= end <= len(sequence):
            seqs.append(sequence[start - 1:end])  # UniProt is 1-indexed
    return seqs if seqs else None

def main(base_dir):
    """
    Main function that processes protein pair directories and generates reference files.
    
    This function:
    1. Loads domain annotation CSV files
    2. Iterates through protein pair directories
    3. Extracts domain sequences for each protein
    4. Creates ref_file_updated.csv for each pair
    
    Args:
        base_dir (str): Base directory containing protein pair subdirectories and annotation CSV files
        
    Note:
        Expects subdirectory names in format 'UNIPROT1_UNIPROT2' and requires
        domain annotation CSV files in the base directory.
    """
    csv1_path = os.path.join(base_dir, "human_interpro_ppi_domains_consolidated_min3.csv")
    csv2_path = os.path.join(base_dir, "human_uniprot_ppi_sites_min3.csv")

    # Load CSVs
    df1 = pd.read_csv(csv1_path)
    df2 = pd.read_csv(csv2_path)

    # Combine domain annotations from both sources
    domain_df = pd.concat([df1[["UniProt_AC", "Start", "End"]], df2[["UniProt_AC", "Start", "End"]]])
    domain_df.dropna(subset=["Start", "End"], inplace=True)
    domain_df["Start"] = domain_df["Start"].astype(int)
    domain_df["End"] = domain_df["End"].astype(int)

    for folder in os.listdir(base_dir):
        folder_path = os.path.join(base_dir, folder)
        if not os.path.isdir(folder_path):
            continue

        if "_" not in folder:  # adjust for your naming format
            continue
 
        a_id, b_id = folder.split("_", 1)
        a_id = a_id.upper()
        b_id = b_id.upper()

        a_seqs = extract_domains(a_id, domain_df)
        b_seqs = extract_domains(b_id, domain_df)

        if not a_seqs or not b_seqs:
            print(f"[SKIP] No domains found for: {a_id if not a_seqs else ''} {b_id if not b_seqs else ''}")
            continue

        if len(a_seqs) < len(b_seqs):
            a_seqs += [a_seqs[-1]] * (len(b_seqs) - len(a_seqs))
        elif len(b_seqs) < len(a_seqs):
            b_seqs += [b_seqs[-1]] * (len(a_seqs) - len(b_seqs))

        out_df = pd.DataFrame({
            "Sequences Chain A": a_seqs,
            "Sequences Chain B": b_seqs
        })

        out_path = os.path.join(folder_path, "ref_file_updated.csv")
        out_df.to_csv(out_path, index=False)
        print(f"[SAVED] {out_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate domain reference files from UniProt sequences.")
    parser.add_argument("--base_dir", required=True, help="Base directory containing protein pair folders and CSV files.")
    args = parser.parse_args()
    main(args.base_dir)
