# pli-analyzer
Protein-ligand interaction analyzer

## Overview

This tool extracts the interaction interface of protein, RNA, and DNA complexes and can compare it to known sequences of binding domains/motifs or to other structures.

## Installation

1. Install the following dependencies: [miniconda3](https://docs.anaconda.com/miniconda/)
2. Clone the repository.
```
git clone https://github.com/ChangLab/pli-analyzer.git
```
3. Navigate to the directory
```
cd pli-analyzer
```

4. Install dependencies

```
conda env create -f environment.yml
```

5. Activate the environment:
```
conda activate pli-analyzer
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
   
## Usage

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


## Epitope Analysis

Extracts the solvent-accessible surface area (SASA), total surface area, and volume for specified epitope sequence ranges.

### Installation and Usage

1. Install ChimeraX (https://www.cgl.ucsf.edu/chimerax/download.html)
2. Open epitope_analysis.py and edit the variable in the first line (input_dir) to be the file path to the input directory
3. In the input directory, create an additional text file named epitope.txt. This text file will contain the sequence ranges you want to query. The first line of the file should contain the sequence chain, for example "/A" or "/B". The following lines should contain each of the ranges, separated by a dash with no spaces (for example, "4-13" or "20-27"). See an example file below.
<img width="216" alt="Screen Shot 2024-09-02 at 8 11 40 PM" src="https://github.com/user-attachments/assets/659cc0b6-d194-42e8-ad10-847e7d870015">
<br>
4. Open ChimeraX
5. In the top menu, select File -> Open -> epitope_analysis.py
6. A CSV file will be created in the input directory, outputting the SASA, total surface area, and volume for the sequence ranges across all input PDB files.
<img width="582" alt="Screen Shot 2024-09-02 at 8 16 40 PM" src="https://github.com/user-attachments/assets/3c1fa7a2-1268-405b-90f2-b27c6885627e">


  









