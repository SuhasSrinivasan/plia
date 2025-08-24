# PLIA
Protein-ligand interaction analyzer

## Overview

This tool extracts the interaction interface of protein, RNA, and DNA complexes and can compare it to known sequences of binding domains/motifs or to other structures.

## Installation

1. Install the following dependencies: [miniconda3](https://docs.anaconda.com/miniconda/)
2. Clone the repository.
```
git clone https://github.com/ChangLab/plia.git
```
3. Navigate to the directory
```
cd plia
```

4. Install dependencies

```
conda env create -f environment.yml
```

5. Activate the environment:
```
conda activate plia
```

6. Install Voronota

### Installing Voronota

Voronota calculates atom contacts through the Voronoi diagram of atomic balls. 

Access the Github Repo here https://github.com/kliment-olechnovic/voronota

#### Quick Install Guide

1. Download the latest package
```
wget https://github.com/kliment-olechnovic/voronota/releases/download/v1.29.4198/voronota_1.29.4198.tar.gz

OR if wget is not installed use the below command:

curl -LO https://github.com/kliment-olechnovic/voronota/releases/download/v1.29.4198/voronota_1.29.4198.tar.gz
```
2. Unpack the package
```
tar -xf voronota_1.29.4198.tar.gz
rm voronota_1.29.4198.tar.gz
```
3. Change to the package directory
```
cd ./voronota_1.29.4198
```
4. Run CMake
```
cmake . -DEXPANSION_JS=ON -DEXPANSION_LT=ON
```
5. Compile everything
```
make
```
6. Install everything
```
sudo make install
```
   
# PLI-Analyzer Usage

## Overview
This wrapper script `master_runner.py` automates the workflow for processing protein complex structures and analyzing their interaction interfaces with PLI-analyzer. It orchestrates several helper scripts in sequence, ensuring that input and output files are handled correctly across subdirectories.  

The pipeline expects a base directory with subdirectories for each protein complex to process, named in the format:

`uniprotid1-uniprotid2`

Each subdirectory is processed independently, and results are consolidated at the base directory level.

## Workflow

1. `cif_to_pdb.py`

Converts `.cif` files to `.pdb`.

2. `new_generate_ref.py`

Generates reference files (`ref_file_updated.csv`) inside each subdirectory, representing known interaction sites for the given interactor. These interaction sites are determined by cross-referencing the Uniprot ID with Uniprot and InterPro annotations, represented by the CSVs `human_interpro_ppi_domains_consolidated_min3.csv` and `human_uniprot_ppi_sites_min3.csv`

3. `extract_interface.py`

For each subdirectory:

- Runs interface extraction with Voronota.
- Filters sequences by interaction length, default = 3.

Example command executed:

python3 orchestrator.py
--input_dir <subdir>
--output_dir <subdir>
--path_to_voronota <voronota_path>
--ref_file <subdir>/ref_file_updated.csv
[--inter_output]



5. `summarize_sequences.py`

Aggregates results into a summary file at the base directory.

python3 summarize_sequences.py <base_dir> --output_csv summary.csv


6. `extra_summary.py`

Performs an additional summarization step for related information on IDR % and length.


### `master_runner.py` Usage

```
python3 orchestrator.py
--base_dir /path/to/base/directory
--voronota_path /path/to/voronota
[--min_interaction_length 3]
[--keep_extract_interface]
[--inter_output]
```


### Required Arguments

- `--base_dir`  
  Path to the base directory containing subdirectories for each complex.

- `--voronota_path`  
  Path to the Voronota executable.

### Optional Arguments

- `--min_interaction_length` (default: `3`)  
  Minimum interaction length. The wrapper internally adds `+2` before passing the value to `filter_sequences.py`.

- `--keep_extract_interface`  
  If set, `extract_interface.py` will remain in each subdirectory after execution.

- `--inter_output`  
  If set, adds `--inter_output` to each `extract_interface.py` run. This generates intermediate CSV files with residue-level distances.

## Outputs

### Per Subdirectory

- `ref_file_updated.csv`
- Extracted interface results
- (Optional) Intermediate residue-level CSVs (`--inter_output`)

### Base Directory

- `summary.csv` (from `summarize_sequences.py`)
- Additional outputs from `extra_summary.py`

## Example

```
python3 orchestrator.py
--base_dir /Users/user/pli-analyzer/output
--voronota_path /Users/user/bin/voronota
--min_interaction_length 5
--inter_output
```


### Notes

- Subdirectories **must follow the naming convention**:
- All helper scripts must be available in the **same directory** as the wrapper (`orchestrator.py`).
- The wrapper stops execution if any required reference file (`ref_file_updated.csv`) is missing in a subdirectory.






### As standalone extract_interface.py script

### Input Parameters

| Parameter    | Description                              |
|--------------|------------------------------------------|
| `input_dir` | Path to the input folder with the PDF/CIF files.                   |
| `output_dir`  | Path to where the output files should be created.       |
| `path_to_voronota`    | Path to the directory where Voronota is installed. |
| `ref_file`    | Name of CSV file containing known/putative interacting sequences (optional, see example file). |
| `area`    | Minimum interaction area threshold for the Voronoi diagram of atomic balls (optional, default = 0.01) |
| `distance`    | Distance threshold (in Angstroms) for interacting residues (optional, default = None). |
| `padding`    | Number of residues added as padding on either side of the actual interacting residue (optional, default = 1). |
| `inter_output`    | Create intermediate CSV files with residue-level distances (optional, default = False). |
   
### Example Usage

```
python extract_interface.py --input_dir <path_to_dir> --output_dir <path_to_dir> --path_to_voronota <path_to_voronota> --ref_file <path_to_ref_file> --padding 2 
```

Note: Voronota currently only supports PDB files— compatibility with CIF files will be added to this code later.
   
### Reference File Format

Each column of the CSV corresponds to a protein chain in the PDB/CIF files.   
<img width="784" alt="Screen Shot 2024-07-11 at 11 47 13 AM" src="https://github.com/user-attachments/assets/caa185f3-4b6d-4521-9d9a-9457033cf79d">
   
   
### Output
   
Intermediate CSV: A file is produced for each input PDB/CIF file, and can be toggled ON by the `inter_output` argument.   
<img width="715" alt="Screen Shot 2024-08-12 at 10 44 24 AM" src="https://github.com/user-attachments/assets/5a6f1be1-44d0-48c3-8217-6882abd2e716">

Final CSV with interacting sequences.      
<img width="638" alt="Screen Shot 2024-08-12 at 10 51 42 AM" src="https://github.com/user-attachments/assets/a34a3110-962e-4953-b5b2-aa2ecaa445b3">


