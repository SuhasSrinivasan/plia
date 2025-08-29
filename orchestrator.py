"""
PLIA Pipeline Orchestrator

This script orchestrates the complete protein-ligand interaction analysis pipeline by running
all necessary scripts in sequence. It handles the full workflow from structure conversion
to final data summarization.

Pipeline Steps:
1. Convert CIF files to PDB format
2. Generate reference domain files
3. Extract protein interfaces using Voronota
4. Filter sequences by length
5. Summarize sequence data
6. Enhance with UniProt information

Usage:
    python orchestrator.py --base_dir /path/to/data --voronota_path /path/to/voronota
                          --min_interaction_length 3

Requirements:
    - All PLIA pipeline scripts in the same directory
    - Voronota software installed
    - Internet connection for UniProt data

Author: PLIA Project
"""

import os
import subprocess
import argparse
import shutil

def run_command(cmd, description):
    """
    Execute a shell command with error handling and status reporting.
    
    Args:
        cmd (str): Shell command to execute
        description (str): Description of the operation for user feedback
        
    Raises:
        SystemExit: If command fails, exits the program with error status
    """
    print(f"🟡 {description}...")
    try:
        subprocess.run(cmd, shell=True, check=True)
        print(f"✅ Done: {description}\n")
    except subprocess.CalledProcessError as e:
        print(f"❌ Failed: {description}")
        print(e)
        exit(1)

def main():
    """
    Main orchestrator function that coordinates the entire pipeline.
    
    This function:
    1. Parses command line arguments
    2. Executes each pipeline step in sequence
    3. Handles file copying and cleanup
    4. Provides progress feedback to user
    
    Command line arguments:
        --base_dir: Root directory containing protein pair subdirectories
        --voronota_path: Path to Voronota executable
        --min_interaction_length: Minimum sequence length for interactions (default: 3)
        --keep_extract_interface: Keep temporary script copies in subdirectories
    """
    parser = argparse.ArgumentParser(description="Run full protein interface analysis pipeline.")
    parser.add_argument('--base_dir', required=True, help="Path to base directory containing subdirectories.")
    parser.add_argument('--voronota_path', required=True, help="Path to Voronota executable.")
    parser.add_argument('--min_interaction_length', type=int, default=3,
                        help="Minimum interaction length (default: 3).")
    parser.add_argument('--keep_extract_interface', action='store_true',
                        help="Keep extract_interface.py in subdirectories after execution.")

    args = parser.parse_args()

    base_dir = os.path.abspath(args.base_dir)
    voronota_path = os.path.abspath(args.voronota_path)
    min_length = args.min_interaction_length + 2

    # Step 1: Run cif_to_pdb.py
    run_command(
        f"python3 cif_to_pdb.py --input_dir \"{base_dir}\"",
        "Running cif_to_pdb.py"
    )

    # Step 2: Run new_generate_ref.py
    run_command(
        f"python3 new_generate_ref.py --base_dir \"{base_dir}\"",
        "Running new_generate_ref.py"
    )

    # Step 3: Run extract_interface.py in each subdirectory
    subdirs = [
        os.path.join(base_dir, d)
        for d in os.listdir(base_dir)
        if os.path.isdir(os.path.join(base_dir, d))
    ]

    print(f"🧩 Found {len(subdirs)} subdirectories to process with extract_interface.py\n")

    for subdir in subdirs:
        ref_file_path = os.path.join(subdir, "ref_file_updated.csv")
        if not os.path.exists(ref_file_path):
            print(f"❌ Error: ref_file_updated.csv not found in {subdir}")
            exit(1)

        target_script = os.path.join(subdir, "extract_interface.py")
        shutil.copy("extract_interface.py", target_script)

        cmd = (
            f"python3 \"{target_script}\" "
            f"--input_dir \"{subdir}\" "
            f"--output_dir \"{subdir}\" "
            f"--path_to_voronota \"{voronota_path}\" "
            f"--ref_file \"{ref_file_path}\""
        )

        run_command(cmd, f"Running extract_interface.py in {os.path.basename(subdir)}")

        if not args.keep_extract_interface:
            os.remove(target_script)
            print(f"🗑️  Deleted extract_interface.py from {subdir}\n")

    # Step 4: Run filter_sequences.py
    run_command(
        f"python3 filter_sequences.py --root_dir \"{base_dir}\" --min_length {min_length}",
        "Running filter_sequences.py"
    )

    # Step 5: Run summarize_sequences.py
    summary_csv_path = os.path.join(base_dir, "summary.csv")
    run_command(
        f"python3 summarize_sequences.py \"{base_dir}\" --output_csv \"{summary_csv_path}\"",
        "Running summarize_sequences.py"
    )

    # Step 6: Run extra_summary.py
    run_command(
        f"python3 extra_summary.py \"{base_dir}\"",
        "Running extra_summary.py"
    )

    print("🎉 All steps completed successfully.")

if __name__ == "__main__":
    main()
