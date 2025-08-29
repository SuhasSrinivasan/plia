"""
Sequence Analysis Summarizer

This script analyzes final_sequences_updated.csv files across subdirectories and generates
comprehensive summary statistics including identity/similarity scores, sequence counts,
and AlphaFold confidence metrics. It combines sequence analysis with structural confidence
data to provide a complete overview of protein interaction predictions.

Usage:
    python summarize_sequences.py /path/to/input/dir --output_csv summary.csv

Requirements:
    - pandas library
    - JSON files with AlphaFold confidence scores (*_summary_confidences_0.json)
    - final_sequences_updated.csv files in subdirectories

Author: PLIA Project
"""

import os
import pandas as pd
import json

def analyze_csv(csv_path, json_path):
    """
    Analyze a single CSV file and extract summary statistics.
    
    Args:
        csv_path (str): Path to final_sequences_updated.csv file
        json_path (str): Path to corresponding AlphaFold confidence JSON file
        
    Returns:
        dict or None: Dictionary containing analysis results with keys:
            - subfolder: Name of the containing subdirectory
            - avg_identity: Average sequence identity score
            - avg_similarity: Average sequence similarity score  
            - num_sequences: Total number of sequences analyzed
            - avg_sequence_length: Average length of sequences
            - pct_identity_X_Y: Percentage of sequences in identity score ranges
            - ipTM: Interaction confidence score from AlphaFold
            - pTM: Overall confidence score from AlphaFold
            
    Note:
        Returns None if analysis fails due to missing data or errors
    """
    try:
        csv_path = str(csv_path)
        df = pd.read_csv(csv_path)
        df = df.dropna(subset=["sequence", "identity", "similarity"])

        avg_identity = df["identity"].mean()
        avg_similarity = df["similarity"].mean()
        num_sequences = len(df)
        avg_seq_len = df["sequence"].apply(len).mean()

        # Calculate identity score distribution percentages
        total = len(df)
        pct_80_100 = ((df["identity"] >= 80) & (df["identity"] <= 100)).sum() / total * 100
        pct_60_80  = ((df["identity"] >= 60) & (df["identity"] < 80)).sum() / total * 100
        pct_40_60  = ((df["identity"] >= 40) & (df["identity"] < 60)).sum() / total * 100
        pct_20_40  = ((df["identity"] >= 20) & (df["identity"] < 40)).sum() / total * 100
        pct_0_20   = ((df["identity"] >= 0)  & (df["identity"] < 20)).sum() / total * 100

        # Extract AlphaFold confidence scores if available
        iptm = ptm = None
        if os.path.exists(json_path):
            with open(json_path, 'r') as f:
                data = json.load(f)
                iptm = data.get("iptm")
                ptm = data.get("ptm")

        return {
            "subfolder": os.path.basename(os.path.dirname(csv_path)),
            "avg_identity": avg_identity,
            "avg_similarity": avg_similarity,
            "num_sequences": num_sequences,
            "avg_sequence_length": avg_seq_len,
            "pct_identity_80_100": pct_80_100,
            "pct_identity_60_80": pct_60_80,
            "pct_identity_40_60": pct_40_60,
            "pct_identity_20_40": pct_20_40,
            "pct_identity_0_20": pct_0_20,
            "ipTM": iptm,
            "pTM": ptm
        }

    except Exception as e:
        print(f"Error processing {csv_path}: {e}")
        return None

def collect_summaries(input_dir):
    """
    Collect and analyze all final_sequences_updated.csv files in subdirectories.
    
    Args:
        input_dir (str): Root directory to search for CSV files
        
    Returns:
        pandas.DataFrame: Combined summary data from all processed subdirectories
                         with columns for identity scores, sequence statistics,
                         and AlphaFold confidence metrics
                         
    Note:
        Searches for files ending with '_summary_confidences_0.json' for AlphaFold data
    """
    summaries = []

    for root, dirs, files in os.walk(input_dir):
        if "final_sequences_updated.csv" in files:
            csv_path = os.path.join(root, "final_sequences_updated.csv")
            # Find the JSON file that ends with _summary_confidences_0.json
            json_path = None
            for f in files:
                if f.endswith("_summary_confidences_0.json"):
                    json_path = os.path.join(root, f)
                    break

            summary = analyze_csv(csv_path, json_path)
            if summary:
                summaries.append(summary)

    return pd.DataFrame(summaries)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Summarize final_sequences.csv + confidence JSON across subfolders."
    )
    parser.add_argument(
        "input_dir", 
        help="Path to directory with subfolders containing final_sequences_updated.csv files"
    )
    parser.add_argument(
        "--output_csv", 
        default="alphafold_summary.csv", 
        help="Output summary CSV file name (default: alphafold_summary.csv)"
    )

    args = parser.parse_args()

    df_summary = collect_summaries(args.input_dir)
    df_summary.to_csv(args.output_csv, index=False)
    print(f"Summary saved to {args.output_csv}")
