import os
import pandas as pd
import json

def analyze_csv(csv_path, json_path):
    try:
        csv_path = str(csv_path)
        df = pd.read_csv(csv_path)
        df = df.dropna(subset=["sequence", "identity", "similarity"])

        avg_identity = df["identity"].mean()
        avg_similarity = df["similarity"].mean()
        num_sequences = len(df)
        avg_seq_len = df["sequence"].apply(len).mean()

        total = len(df)
        pct_80_100 = ((df["identity"] >= 80) & (df["identity"] <= 100)).sum() / total * 100
        pct_60_80  = ((df["identity"] >= 60) & (df["identity"] < 80)).sum() / total * 100
        pct_40_60  = ((df["identity"] >= 40) & (df["identity"] < 60)).sum() / total * 100
        pct_20_40  = ((df["identity"] >= 20) & (df["identity"] < 40)).sum() / total * 100
        pct_0_20   = ((df["identity"] >= 0)  & (df["identity"] < 20)).sum() / total * 100

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


            # json_path = None
            # for f in files:
            #     if f.startswith("confidence_") and f.endswith("_model_0.json"):
            #         json_path = os.path.join(root, f)
            #         break
                

            # print ("json path", json_path)
            summary = analyze_csv(csv_path, json_path)
            if summary:
                summaries.append(summary)

    return pd.DataFrame(summaries)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Summarize final_sequences.csv + confidence JSON across subfolders.")
    parser.add_argument("input_dir", help="Path to directory with subfolders")
    parser.add_argument("--output_csv", default="alphafold_summary.csv", help="Output summary CSV file name")

    args = parser.parse_args()

    df_summary = collect_summaries(args.input_dir)
    df_summary.to_csv(args.output_csv, index=False)
    print(f"Summary saved to {args.output_csv}")
