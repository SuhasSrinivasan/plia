"""
CIF to PDB Converter

This script converts Crystallographic Information Files (CIF) to Protein Data Bank (PDB) format
using Biopython. It can process single files or entire directory trees, automatically converting
all .cif files found and optionally removing the original files after successful conversion.

Usage:
    python cif_to_pdb.py --input_dir /path/to/directory

Requirements:
    - Biopython library
    - Valid CIF files for conversion

Author: PLIA Project
"""

import os
import argparse
from Bio.PDB import MMCIFParser, PDBIO

def convert_cif_to_pdb(cif_path, pdb_path):
    """
    Convert a CIF file to a PDB file using Biopython.
    
    Args:
        cif_path (str): Path to input CIF file
        pdb_path (str): Path where output PDB file will be saved
        
    Raises:
        Exception: If CIF parsing or PDB writing fails
    """
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("structure", cif_path)
    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_path)

def process_directory(base_dir):
    """
    Traverse directory and convert all .cif files to .pdb, then delete the .cif files.
    
    This function recursively walks through all subdirectories starting from base_dir,
    finds all .cif files, converts them to .pdb format, and removes the original
    .cif files upon successful conversion.
    
    Args:
        base_dir (str): Root directory to start searching for CIF files
        
    Note:
        Original CIF files are deleted only after successful conversion
    """
    for root, dirs, files in os.walk(base_dir):
        for file in files:
            if file.lower().endswith(".cif"):
                cif_path = os.path.join(root, file)
                pdb_path = os.path.splitext(cif_path)[0] + ".pdb"

                try:
                    print(f"Converting: {cif_path} → {pdb_path}")
                    convert_cif_to_pdb(cif_path, pdb_path)
                    os.remove(cif_path)
                    print(f"Deleted original CIF file: {cif_path}")
                except Exception as e:
                    print(f"[ERROR] Failed to convert {cif_path}: {e}")

def main():
    """
    Main function that parses command line arguments and processes the directory.
    
    Command line arguments:
        --input_dir (str): Base directory to search for CIF files (required)
    """
    parser = argparse.ArgumentParser(description="Convert all .cif files to .pdb in a directory tree.")
    parser.add_argument(
        "--input_dir",
        required=True,
        help="Base directory to search for CIF files"
    )

    args = parser.parse_args()
    base_dir = args.input_dir

    if os.path.isdir(base_dir):
        process_directory(base_dir)
    else:
        print(f"[ERROR] Invalid directory: {base_dir}")

if __name__ == "__main__":
    main()
