"""
Sequence Length Filter for Final Sequences

This script filters and processes final_sequences.csv files in subdirectories, removing
sequences that don't meet minimum length requirements and applying formatting rules.
It specifically targets model_0 files and removes padding characters from sequences.

Usage:
    python filter_sequences.py --root_dir /path/to/directory --min_length 5

Requirements:
    - pandas library
    - Existing final_sequences.csv files in subdirectories

Author: PLIA Project
"""

import os
import pandas as pd
import argparse

def main(root_dir, min_length):
    """
    Main function that processes final_sequences.csv files in subdirectories.
    
    This function:
    1. Recursively searches for final_sequences.csv files
    2. Filters sequences by minimum length and model type
    3. Removes padding characters (first and last character)
    4. Saves filtered results to final_sequences_updated.csv
    
    Args:
        root_dir (str): Root directory to search for final_sequences.csv files
        min_length (int): Minimum sequence length to retain after filtering
        
    Note:
        Only processes sequences from 'model_0' files and removes first/last characters
        as padding before applying length filter.
    """
    # iterate through subdirectories
    for subdir, _, files in os.walk(root_dir):
        if 'final_sequences.csv' in files:
            csv_path = os.path.join(subdir, 'final_sequences.csv')
            updated_csv_path = os.path.join(subdir, 'final_sequences_updated.csv')
            
            try:
                df = pd.read_csv(csv_path)
                
                if 'sequence' in df.columns and 'file' in df.columns:
                    df['sequence'] = df['sequence'].astype(str)
                    df['file'] = df['file'].astype(str)

                    # Filter rows where sequence length >= min_length and only the model_0
                    filtered_df = df[
                        (df['sequence'].str.len() >= min_length) &
                        (df['file'].str.endswith('model_0'))
                    ].copy()
                    
                    # remove the first and last character of each sequence (padding)
                    filtered_df['sequence'] = filtered_df['sequence'].apply(lambda x: x[1:-1])
                    
                    filtered_df.to_csv(updated_csv_path, index=False)
                    print(f"Updated file written to: {updated_csv_path}")
                else:
                    print(f"Required columns not found in: {csv_path}")
            except Exception as e:
                print(f"Failed to process {csv_path}: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter and update final_sequences.csv files.")
    parser.add_argument("--root_dir", required=True, help="Root directory containing subfolders with final_sequences.csv files.")
    parser.add_argument("--min_length", type=int, required=True, help="Minimum sequence length to keep.")
    args = parser.parse_args()
    
    main(args.root_dir, args.min_length)