import os
import pandas as pd
import argparse

def main(root_dir, min_length):
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