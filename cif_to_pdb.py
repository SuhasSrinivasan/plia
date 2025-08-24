import os
import argparse
from Bio.PDB import MMCIFParser, PDBIO

def convert_cif_to_pdb(cif_path, pdb_path):
    """Convert a CIF file to a PDB file using Biopython."""
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("structure", cif_path)
    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_path)

def process_directory(base_dir):
    """Traverse directory and convert all .cif files to .pdb, then delete the .cif files."""
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
