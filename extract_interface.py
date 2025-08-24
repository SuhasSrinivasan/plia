import Bio.PDB #Imports
from Bio.PDB import *
import argparse
import os
import glob
import numpy as np
import string
import pandas as pd
from Bio.PDB.Polypeptide import protein_letters_3to1
from Bio import BiopythonWarning
import warnings
from Bio.PDB import PDBParser
from Bio.PDB import MMCIFParser, PDBIO, PDBParser
from Bio import Align
import csv
from collections import defaultdict
import re
import subprocess



def is_nucleic_acid(residue):
    """Checks if a residue is a nucleic acid

    Args: 
        residue (residue class object): Accessed through pipeline of structure --> model --> chain --> residue, consists of atoms

    Returns:
        boolean: Whether the residue is a nucleic acid

    note: A, T, G, etc represent ribonucleic acids, while DA, DT, DG, etc represent deoxyribonucleic acids
    """
    nucleic_acids = {"A", "T", "G", "C", "U", "DA", "DT", "DG", "DC", "DU"} #Dictionary of nucleic acids
    return residue.get_resname() in nucleic_acids


def is_amino_acid(residue):
    """Checks if residue is amino acid

    Args:
        residue (residue class): See is_nucleic_acid

    Returns:
        boolean: Whether the residue is an amino acid
    """
    return (residue.get_resname() in standard_aa_names) #standard_aa_names is a dictionary with all amino acid codons (3 letters)




def run_padding(chain, position, pad_distance): #Default 1 on each side
    """Returns a padded sequence of residues around a residue

    Args:
        chain (list of residue objects): Chain A, B, C, etc
        position (int): Position of "central" residue that the padding is built around on either side 
        pad_distance (int): default = 1

    Returns:
        string: Padded Sequence

    
    Example:

    Chain A position 5 (arginine) is passed through, with a pad distance of 3

    Returned sequence: "KGGALPS" --> A is in the center, with the 3 residues around it
    """

    # print("chain[0]", chain[0])
    # print("position", position)

    # Your single residue string
    residue_string = (str(chain[0]))

    pattern = r"resseq=(\d+)"

    match = re.search(pattern, residue_string)

    resseq = match.group(1)
    # print(resseq)

    difference = 0

    if resseq == 1: #check if the residues in the chain are ordered properly: if not, calculate how far off the residues are
        if (position - pad_distance  + 1 < 1): #If position is at the start of the chain
            residues = [chain[res_id] for res_id in range(0, position + pad_distance)]
        elif (position + pad_distance >= len(chain)): #If position is at the end of the chain
            residues = [chain[res_id] for res_id in range(position - pad_distance, len(chain))]
        else:
            residues = [chain[res_id] for res_id in range(position - pad_distance , position + pad_distance + 1)] 
    else:
        updated_pos = position - int(resseq) + 1 #updated, "correct" position
        # difference == int(resseq) - 1
        if (updated_pos - pad_distance  < 0): #If position is at the start of the chain
            residues = [chain[res_id] for res_id in range(0, updated_pos + pad_distance)]
        elif (updated_pos + pad_distance  >= len(chain)): #If position is at the end of the chain
            residues = [chain[res_id] for res_id in range(updated_pos - pad_distance, len(chain))]
        else:
            residues = [chain[res_id] for res_id in range(updated_pos - pad_distance , updated_pos + pad_distance + 1)] 
    

    # Convert residues to one-letter codes and concatenate them into a string
    if len(residues[0].get_resname()) == 3:
        sequence = ''.join([protein_letters_3to1[(residue.get_resname())] for residue in residues])
    elif len(residues[0].get_resname()) == 2: #Checking for if it's DNA, where bases are defined as DA, DC, DG, etc. then, .get_resname can be used instead of having to convert text formats
        sequence = ''.join([(residue.get_resname())[1] for residue in residues]) #Using only the second character, ie C, G of DC, DG
    else: #Only one character, means it's an RNA
        sequence = ''.join([(residue.get_resname()) for residue in residues])
    # print(f"Position {position - pad_distance} to {position + pad_distance + 1} is {sequence}")
    return sequence




def calc_dist_matrix(ch1, ch2, model, distance_threshold, pad_distance, path_to_file, path_to_voronota, area_threshold):
    """Runs VoroContacts to calculate contacts above a certain area threshold for 2 chains in a PDB File

    Args:
        ch1 (chain class object): First chain
        ch2 (chain class object): Second chain
        model (model class object): Model from PDB File
        distance_threshold (float): _description_
        pad_distance (int): Distance for padding on either side, default = 1
        path_to_file (string): Path to PDB file that is being run
        path_to_voronota (string): path to the directory where voronota is installed
        area_threshold (float): minimum area value to consider residues as in contact

    Returns:
        _type_: _description_
        Iterates between all pairs of residues of these two chains and returns various lists of the residue pairs that have distances below the distance threshold


    Returns:
        first residues (list, str): List of the first residues that had a distance with another residue that was less than the distance threshold, each residue is formatted as its number then amino acid, for example 25GLY
        second residues (list, str): List of the corresponding second residues that had a distance with another residue that was less than the distance threshold
        distances (list, int): List of the distances between the residue pairs
        first_pads (list, str): Padded residue sequences for  the first residue (residue one) in a pair
        second_pads (list, str): Padded residue sequences for  the second residue (residue two) in a pair
        first_residues_num (list, str): Only the numerical position of the same residues in first_residues
        second_residues_num (list, str): Only the numerical position of the same residues in second_residues
        area (list, int): List of the area overlap (already filtered to be > area_threshold) between the residue pairs

    """
    print("processing chains", ch1, ch2)
    print("path to file", path_to_file)


    generate_balls = f"./voronota get-balls-from-atoms-file --annotated < {path_to_file} > annotated_balls.txt           "
    calculate_contacts = f"./voronota calculate-contacts --annotated <{path_to_voronota}/annotated_balls.txt> annotated_contacts.txt"
    if distance_threshold != None:
        query_contacts = f"./voronota query-contacts --match-min-area {area_threshold} --match-max-dist {distance_threshold} --no-solvent --no-same-chain <{path_to_voronota}/annotated_contacts.txt> queried_contacts.txt"
    else:
        query_contacts = f"./voronota query-contacts --match-min-area {area_threshold} --no-solvent --no-same-chain <{path_to_voronota}/annotated_contacts.txt> queried_contacts.txt"


    output_file_path = path_to_voronota

    os.chdir(path_to_voronota)

    # Execute the command
    result = subprocess.run(generate_balls, shell=True, capture_output=False, text=True)
    result = subprocess.run(calculate_contacts, shell=True, capture_output=False, text=True)
    result = subprocess.run(query_contacts, shell=True, capture_output=False, text=True)


    unique_pairs = set()

    # Regex pattern
    pattern = r"c<([^>]+)>r<([^>]+)>a<([^>]+)>R<([^>]+)>A<([^>]+)>"

    first_residues = []
    chain_one = []
    chain_two = []
    second_residues = []
    first_residues_num = []
    second_residues_num = []
    area = []
    distances = []
    first_pads = []
    second_pads = []


    with open('queried_contacts.txt', 'r') as file: #open the final outputted file (containing the filtered contacts) from Vorocontacts
        for line in file: #read in every line from the file, where each line represents a pair of contacted residues
            line = line.strip()
            if line:
                parts = line.split()
                match1 = re.match(pattern, parts[0])
                match2 = re.match(pattern, parts[1])

                if match1 and match2:
                    values1 = match1.groups()
                    values2 = match2.groups()

                    values1 = list(values1)
                    values2 = list(values2)

                    del values1[2]
                    del values1[-1]

                    del values2[2]
                    del values2[-1]

                    # print("Values from string1:", values1)
                    # print("Values from string2:", values2)

                    res_pair = []
                    res_pair.extend(values1)
                    res_pair.extend(values2)

                    # print("final pair", res_pair)

                    #check if chains match query
                    if ((ch1 == values1[0]) and (ch2 == values2[0])):
                        if tuple(res_pair) not in unique_pairs:
                            # final_pairs.append(res_pair)
                            chain_one.append(res_pair[0])
                            first_residues_num.append(int(res_pair[1]) - 1)
                            first_residues.append(f"{res_pair[1]}{res_pair[2]}")
                            chain_two.append(res_pair[3])
                            second_residues_num.append(int(res_pair[4]) - 1)
                            second_residues.append(f"{res_pair[4]}{res_pair[5]}")
                            area.append(parts[2])
                            distances.append(parts[3])
                            # print("model chain", model[res_pair[0]])

                            chain1_res = []
                            nucleic_acids = {"A", "T", "G", "C", "U", "DA", "DT", "DG", "DC", "DU"}

                            for residue in model[res_pair[0]]:
                                if (residue.get_resname() in standard_aa_names) or (residue.get_resname() in nucleic_acids):
                                    chain1_res.append(residue)
                            

                            chain2_res = []
                            for residue in model[res_pair[3]]:
                                if (residue.get_resname() in standard_aa_names) or (residue.get_resname() in nucleic_acids):
                                    chain2_res.append(residue)

                            # print("chain one res", chain1_res)
                            # print("len of ch1_res", len(chain1_res))
                            # print("chain two res", chain2_res)
                            # print("line queried, ", line)
                            first_pads.append(run_padding(chain1_res, int(res_pair[1]) - 1, pad_distance))
                            second_pads.append(run_padding(chain2_res, int(res_pair[4]) - 1, pad_distance))

                            unique_pairs.add(tuple(res_pair))


    # # first_residues, second_residues, answer, distances, first_pads, second_pads, first_residues_num, second_residues_num
    # print(f"chain_one: {chain_one}")
    # print(f"first_residues_num: {first_residues_num}")
    # print(f"first_residues: {first_residues}")
    # print(f"chain_two: {chain_two}")
    # print(f"second_residues_num: {second_residues_num}")
    # print(f"second_residues: {second_residues}")
    # print(f"area: {area}")
    # print(f"distances: {distances}")


    return first_residues, second_residues,  distances, first_pads, second_pads, first_residues_num, second_residues_num, area






def create_final_columns(first_residues_num, first_pads):
    """Synthesizes all of the padded ranges into final "binding interfaces"

    Args:
        first_residues_num (list): List of residue indicies 
        first_pads (list): List with all the padded sequences

    Returns:
        final_ranges (list): List containing "final strings" that represent the synthesized ranges

    Note: The final ranges are synthesized by combining any of the padded residue ranges that are continuous. For example, if the residues 8, 9, and 10 
    have a residue pair < the distance threshold with sequences ANNL, NNLG, and NLGG, the final range would be ANNLGG

    """
    for i, res in enumerate(first_residues_num):
        first_residues_num[i] = int(res)

    final_ranges = []
    current_string = ""
    for i, res in enumerate(first_residues_num): 
        if i == 0:
            num_consec = 1 #Keeping track of the # of consecutive residues
            current = first_residues_num[0]
            current_string = first_pads[0]


        if i >= len(first_residues_num) - 1: #If we have reached the end of the list of residues 

            final_ranges.append(current_string) #Append whatever current string there is so far       
            return final_ranges


        if first_residues_num[i] == first_residues_num[i + 1]: #If the next residue is equal to the current residue (for example a residue might have <10 A with 2 different residues on the other chain)
            num_consec += 1
        elif first_residues_num[i] == first_residues_num[i+1] - 1: #Directly above
            if len(first_pads[i]) <= len(first_pads[i+1]):
                num_consec = 1
                current += 1
                current_string = current_string + (first_pads[i+1])[-1] 
            else:
                num_consec = 1
                current += 1
        else: ## Big jump for example 1 --> 10, where it should be considered as seperate final ranges

            num_consec = 1 #Resetting
            current = first_residues_num[i+1]
            final_ranges.append(current_string)
            # print("Final string", current_string)
            current_string = first_pads[i+1]



def calculate_similarity(score, seq2):
    """Given the "score" of the best alignment, returns the percentage of bases that are an exact match (what the score is score)


    Args:
        score (int): Score of the best alignment generated by the aligner -> alignments[0]
        seq2 (string): Shorter sequence that is being aligned to the binding site

    Returns:
        float: Percentage sequence similarity
    """

    length = (len(seq2)) #calculates based off of the 
    return (score / length) * 100



def sequence_alignments(binding_sites, seq2): #seq1 is the longer sequence, seq2 is the shorter sequence
    """Given the binding site and comparison sequence, return the alignment with the highest % sequence similarity. note: gaps are allowed and not penalized

    Args:
        binding_sites (list <string>): Known binding sites from the input csv
        seq2 (string): Shorter sequence that is being aligned to the binding site

    Returns:
        float: Highest similarity score (as percentage, so 0 - 100)
    """

    scores = []
    aligner = Align.PairwiseAligner() #Creating the aligner from Biopython

    aligner.mode = 'global'
    aligner.match_score = 1 #1 point for exact
    aligner.mismatch_score = 0 #No penalty for mismatching bases, no gap penalty


    for seq1 in binding_sites: #Iterating through all of the binding sites

        if seq1 == "":
            scores.append(0)
        else:
            alignments = aligner.align(seq1, seq2)

        # Print the alignment in a nice format
        # for alignment in alignments:
        #     print(format_alignment(*alignment))

        # Calculate similarity for the best alignment
            best_alignment = alignments[0]
            similarity_percentage = calculate_similarity(best_alignment.score, seq2)
            scores.append(similarity_percentage)

            # print(f"Percentage sequence similarity: {similarity_percentage:.2f}%")

    return max(scores)




def identity_score(binding_sites, seq2):
    """Calculates the identity score between the input sequence and inputted binding sites (different from similarity scores in that no gaps are allowed)

    Args:
        binding_sites (list <string>): Binding sites from the input csv
        seq2 (string): Sequence that is being compared to the binding site

    Returns:
        float: Percentage identity score (0 - 100)
    """
    scores = []
    for k, seq1 in enumerate(binding_sites):
        if len(seq1) < len(seq2):
            seq1, seq2 = seq2, seq1 #switch so that seq1 is the longer, "global" ref seq
        
        max_score = float('-inf')
        best_alignment = None

        # Slide the shorter sequence along the longer sequence, sort of like a frameshift
        for i in range(len(seq1) - len(seq2) + 1):
            score = 0
            for j in range(len(seq2)):
                if seq1[i + j] == seq2[j]:
                    score += 1  # Each match contributes a score of 1

            if score > max_score:
                max_score = score
                best_alignment = i  

            
        scores.append(max_score)
            
        # print("Max score", max_score)
    # print("Scores", scores)
    best_site = scores.index(max(scores))

    if (max(scores) / len(seq2)) * 100 > 100:
        return 100, binding_sites[best_site]
    else:
        return (max(scores) / len(seq2)) * 100, binding_sites[best_site]



def aa_residues(chain):

    """Removes residues that aren't a nucleic or amino acid

    Args:
        chain (chain class object): Chain A, B, etc

    Returns:
        (list): List of only the amino acid + nucleic acid residues
    """
    # print("am i running these??")

    nucleic_acids = {"A", "T", "G", "C", "U", "DA", "DT", "DG", "DC", "DU"}

    aa_only = []
    for i in chain:
        if (i.get_resname() in standard_aa_names) or (i.get_resname() in nucleic_acids):
            # print(i, "getting added")
            aa_only.append(i)
    return aa_only




#####################

def extract_structure(file_path, file_type):
    """Given a file path and type, it uses the corresponding parser to extract the structure object

    Args:
        file_path (string): Path to PDB/CIF file
        file_type (string): CIF or PDB

    Raises:
        ValueError: file type besides CIF or PDB

    Returns:
        Structure class object: Biopython object in the hierarchy of Structure --> Model  -> Chain --> Residue --> Atom
    """
    if file_type == 'cif':
        parser = MMCIFParser(QUIET=True)
    elif file_type == 'pdb':
        parser = PDBParser(QUIET=True)
    else:
        raise ValueError("Unsupported file type: only '.cif' and '.pdb' are supported.")
    
    structure_id = os.path.splitext(os.path.basename(file_path))[0]
    structure = parser.get_structure(structure_id, file_path)
    return structure

def remove_empty_values(dict):
    """removes empty list values in a dictionary

    Args:
        dict (dictionary): dictionary of lists

    Returns:
        dictionary: dictionary with blank values removed
    """
    return {k: [v for v in val if v != ""] for k, val in dict.items()}


def create_csv(model, chain_number_1, chain_number_2, pdb_name, file_path, distance_threshold, pad_distance, binding_sites, path_to_voro, area_threshold):
    """Puts all of the elements together to create a csv between two chains

    Args:
        model (model class): Subunit of a structure
        chain_number_1 (chain class): Chain 1
        chain_number_2 (chain class): Chain 2
        pdb_name (string): PDB name, for example "5tph"
        distance_threshold (int): Threshold for residue distances
        pad_distance (int): Distance for padding on both sides
        binding_sites (list <strings>): Binding sites from input csv

    Returns:
        df (Pandas dataframe): Intermediate CSV for each PDB/CIF file
        df_2 (Pandas dataframe): Sequence range information that is synthesized into the final_sequences.csv
    """

    d = dict(enumerate(string.ascii_uppercase, 1))
    ch1 = (d[chain_number_1]) #Converts the number to string, for example 1 -> A, 2 -> B, 3 -> C, etc
    ch2 = (d[chain_number_2])


    chain_a = aa_residues(model[ch1]) 
    chain_b = aa_residues(model[ch2])

    first_residues, second_residues, distances, first_pads, second_pads, first_residues_num, second_residues_num, area = calc_dist_matrix(ch1, ch2, model, distance_threshold, pad_distance, file_path, path_to_voro, area_threshold)

    if binding_sites != None:
        binding_sites = remove_empty_values(binding_sites) #Account for unequal or empty column lengths for inputted binding sites

    for i in range (len(second_residues_num)):
        second_residues_num[i] = int(second_residues_num[i])


    zipped_lists = list(zip(second_residues_num, second_pads))

    sorted_zipped_lists = sorted(zipped_lists, key = lambda x: x[0])

    res = [[i for i, j in sorted_zipped_lists],
       [j for i, j in sorted_zipped_lists]]

    second_residues_num = list(res[0])
    second_pads_ = list(res[1])

    total_length = len(first_pads)

    final_ranges_one = (create_final_columns(first_residues_num, first_pads))

    sim_scores_1 =  []
    id_scores_1 = []
    bind_sites_1 = []
    for bind_range in final_ranges_one:
        if binding_sites != None:
            sim_scores_1.append(round(sequence_alignments(binding_sites[f"Sequences Chain {ch1}"], bind_range), 2))
            id_num, bind_site = identity_score(binding_sites[f"Sequences Chain {ch1}"], bind_range)
            id_scores_1.append(round(id_num, 2))
            bind_sites_1.append(bind_site)

    sim_scores_1_copy = sim_scores_1[:]
    id_scores_1_copy = id_scores_1[:]
    final_ranges_1_copy = final_ranges_one[:]


    if final_ranges_one != None:
        if len(final_ranges_one) < total_length:
            for i in range (total_length - len(final_ranges_one)):
                final_ranges_one.append("")
                sim_scores_1.append("")
                id_scores_1.append("")


    #Final ranges for second residue in pairs

    final_ranges_two = (create_final_columns(second_residues_num, second_pads_))


    sim_scores_2 = []
    id_scores_2 = []
    bind_sites_2 = []

    for bind_range in final_ranges_two:
        if binding_sites != None:
            print("binding sites", binding_sites)
            sim_scores_2.append(round(sequence_alignments(binding_sites[f"Sequences Chain {ch2}"], bind_range), 2))
            id_num, bind_site = identity_score(binding_sites[f"Sequences Chain {ch2}"], bind_range)
            id_scores_2.append(round(id_num, 2))
            bind_sites_2.append(bind_site)

    sim_scores_2_copy = sim_scores_2[:]
    id_scores_2_copy = id_scores_2[:]
    final_ranges_2_copy = final_ranges_two[:]

    if final_ranges_two != None:
        if len(final_ranges_two) < total_length:
            for i in range (total_length - len(final_ranges_two)):
                final_ranges_two.append("") #Pandas only accepts arrays of equal length in the dataframe, so add blank elements to reach equal length
                sim_scores_2.append("")
                id_scores_2.append("")
        
    chain_one = [f"{ch1}"] * len(first_residues)
    chain_two= [f"{ch2}"] * len(first_residues)

    # import pylab
    # pylab.matshow(numpy.transpose(dist_matrix))
    # pylab.colorbar()
    # pylab.show()


    max_length = max(len(first_residues), len(second_residues))
    pdb_names = [f"{pdb_name}"] * max_length


    if binding_sites != None:
        main_data = { #First data frame (one for each pdb/cif input filie)
            "file": pdb_names,
            "chain_1": chain_one,
            "chain_2": chain_two,
            "residue_1": first_residues,
            "residue_2": second_residues,
            "area": area,
            "distance": distances,
            "pads_1": first_pads,
            "pads_2": second_pads,
            "sequence_1": final_ranges_one,
            "sequence_2": final_ranges_two,

        }

        data_scores = { #Second data file that synthesizes all the binding ranges across all of the input files 
        "Final ranges for residue one": final_ranges_one,
        "Similarity scores for res 1": sim_scores_1,
        "Identity scores res1": id_scores_1,
        "Final ranges for residue two": final_ranges_two,
        "Similarity scores for res 2": sim_scores_2,
        "Identity scores res2": id_scores_2
        }

        df = pd.DataFrame(main_data)
        df_2 = pd.DataFrame(data_scores)


        print("Data frame between chains", ch1, ch2)
        print (df)


        return df, df_2, sim_scores_1_copy , id_scores_1_copy, final_ranges_1_copy, sim_scores_2_copy , id_scores_2_copy, final_ranges_2_copy, bind_sites_1, bind_sites_2

    else:

        main_data = { #First data frame (one for each PDB/CIF input filie)
            "file": pdb_names,
            "chain_1": chain_one,
            "chain_2": chain_two,
            "residue_1": first_residues,
            "residue_2": second_residues,
            "area": area,
            "distance": distances,
            "pads_1": first_pads,
            "pads_2": second_pads,
            "sequence_1": final_ranges_one,
            "sequence_2": final_ranges_two,
        }




        df = pd.DataFrame(main_data)



        print("Data frame between chains", ch1, ch2)
        print (df)


        return df, None, sim_scores_1_copy , id_scores_1_copy, final_ranges_1_copy, sim_scores_2_copy , id_scores_2_copy, final_ranges_2_copy, bind_sites_1, bind_sites_2





def main(input_dir, output_dir, path_to_voro, binding_file_name, area_threshold, distance_threshold, pad_distance, inter_output):
    """Iterates through all appropriate input files and runs them, synthesizing sequence ranges together to make the final sequence CSV  

    Args:
        input_dir (string): File path to input directory
        output_dir (string): File path to output directory (Will make the folder if it doesn't exist)
        path_to_voro (string): Path to directory where Voronota is installed
        binding_file_name (string): Path to CSV file containing known interaction sequences
        area_threshold (float): area threshold for overlap between Voronoi tesselations to consider two residues as in contact.
        distance_threshold (float): Distance threshold (in Angstroms) for interacting residues (optional, default = None).
        pad_distance (int): Number of residues added as padding on either side of the actual interacting residue (optional, default = 1).
        inter_output (bool): Does not create intermediate CSV files (optional, default = False).
    """


    # Ensure input directory exists
    if not os.path.isdir(input_dir):
        print(f"Error: The input directory '{input_dir}' does not exist.")
        return
    
    # Create the output directory if it does not exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Find all .cif and .pdb files in the input directory
    cif_files = glob.glob(os.path.join(input_dir, '*.cif'))
    pdb_files = glob.glob(os.path.join(input_dir, '*.pdb'))
    
    all_files = cif_files + pdb_files

    columns = defaultdict(list) # each value in each column is appended to a list

    

    # print(os.getcwd())
    if binding_file_name != None:
        with open(binding_file_name) as f:
            reader = csv.DictReader(f) # Read rows into a dictionary format
            for row in reader: # Read a row as {column1: value1, column2: value2,...}
                for (k,v) in row.items(): # Go over each column name and value 
                    columns[k].append(v) # Append the value into the appropriate list
    


        binding_sites = columns
    else:
        columns = None
    final_sequence_scores = []

    if not all_files:
        print(f"No .cif or .pdb files found in the directory '{input_dir}'.")
        return
    
    ranges_1 = []
    sim_scores_1= []
    id_scores_1 = []
    ranges_2 = []
    sim_scores_2 = []
    id_scores_2 = []
    file_names = []
    chain_names = []
    bind_sites = []


    all_lists = []

    d = dict(enumerate(string.ascii_uppercase, 1)) #Dictionary to convert numbers to letters, 1 -> A, 2 -> B

    #Process each file
    for file_path in all_files:
        file_type = os.path.splitext(file_path)[1][1:]  # Get the file extension without the dot
        structure = extract_structure(file_path, file_type)
        pdb_name = os.path.splitext(os.path.basename(file_path))[0]
        data_frames = []

        model = structure[0]
        num_chains = len(list(model))

        print("Processing file:", file_path)
        print("Number of chains:", num_chains)
        for i in range(num_chains):
            for j in range(i + 1, num_chains):
                d = dict(enumerate(string.ascii_uppercase, 1))
                ch1 = (d[i + 1]) #Converts the number to string, for example 1 -> A, 2 -> B, 3 -> C, etc
                ch2 = (d[j + 1])
                if aa_residues(model[ch1]) != [] and aa_residues(model[ch2]):
                    current_df, df_2, sim_scores_1_copy , id_scores_1_copy, final_ranges_1_copy, sim_scores_2_copy , id_scores_2_copy, final_ranges_2_copy, bind_sites_1, bind_sites_2  = create_csv(structure[0], i + 1, j + 1, pdb_name, file_path, distance_threshold, pad_distance, columns, path_to_voro, area_threshold)
                    data_frames.append(current_df)
                    ranges_1.extend(final_ranges_1_copy)
                    ranges_1.extend(final_ranges_2_copy)
                    sim_scores_1.extend(sim_scores_1_copy)
                    sim_scores_1.extend(sim_scores_2_copy)
                    id_scores_1.extend(id_scores_1_copy)
                    id_scores_1.extend(id_scores_2_copy)
                    if binding_file_name != None:
                        bind_sites.extend(bind_sites_1)
                        bind_sites.extend(bind_sites_2)

                    file_names.extend([f"{pdb_name}"] * (len(final_ranges_1_copy) + len(final_ranges_2_copy)))
                    chain_names.extend([d[i + 1]]*len(final_ranges_1_copy)) 
                    chain_names.extend([d[j + 1]]*len(final_ranges_2_copy)) 

        if (inter_output == True):
            os.chdir(output_dir)
            print("data frames", data_frames)
            result = pd.concat(data_frames)
            result.to_csv(f"{pdb_name}_binding_interface.csv",  index=False)

    if binding_file_name != None:


        all_lists.append(file_names)
        all_lists.append(chain_names)
        all_lists.append(ranges_1)
        all_lists.append(id_scores_1)
        all_lists.append(sim_scores_1)
        all_lists.append(bind_sites)

        os.chdir(output_dir)

        max_length = max(len(lst) for lst in all_lists)

        #Pad lists with empty strings
        padded_lists = [lst + [""] * (max_length - len(lst)) for lst in all_lists]

        headers = ["file", "chain", "sequence", "identity", "similarity", "corresponding reference site"] 

        sequence_data = pd.DataFrame({f'col_{i+1}': padded_lists[i] for i in range(len(padded_lists))})

        sequence_data.columns = headers

        sequence_data = sequence_data.sort_values(by="identity", ascending=False)

        sequence_data.to_csv(f"final_sequences.csv",  index=False)




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process .pdb and .cif files in a directory.')
    
    parser.add_argument('--input_dir', type=str, required=True, help='The input file directory containing .pdb files')
    parser.add_argument('--output_dir', type=str, required=True, help='The output file directory')
    parser.add_argument('--path_to_voronota', type=str, required=True, help='The file directory where voronota is installed')
    parser.add_argument('--ref_file', type = str, required = False, help = 'CSV with list of binding sites (sequences)')
    parser.add_argument('--area', type=float, default = 0.01, required=False, help='minimum area threshold for VDW interactions')
    parser.add_argument('--distance', type=float, default = 10, required=False, help='distance threshold in angstroms')
    parser.add_argument('--padding', type=int, default = 1, required=False, help='# of residues added as "padding" on either side of the residues')
    parser.add_argument('--inter_output', type=bool, default = False, required=False, help='creates intermediate csvs: will create by default')


    args = parser.parse_args()
    
    main(args.input_dir, args.output_dir, args.path_to_voronota, args.ref_file, args.area, args.distance, args.padding, args.inter_output)



