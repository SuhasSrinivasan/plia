input_dir = "/Users/fralian/Desktop/test"

import os
import chimerax
from chimerax.core.commands import run as rc  # use 'rc' as shorthand for runCommand
from io import StringIO

os.chdir(input_dir)

file_names = [fn for fn in os.listdir(".") if fn.endswith(".pdb")]

output_file_name = "epitope_results.csv"

with open(output_file_name, "w") as output_file:

    for fn in file_names:
        epitope_file_name = "epitope.txt"
        ranges = []

        # Read the chain format and residue ranges from the epitope file
        with open(epitope_file_name, 'r') as file:
            chain_format = file.readline().strip()
            print("Chain format:", chain_format)
        
            for line in file:
                range_start, range_end = map(int, line.strip().split('-'))
                ranges.append((range_start, range_end))
                
        print("Ranges:", ranges)

        # Open the PDB file
        rc(session, "open " + fn)

        # output_file_name = f"{os.path.splitext(fn)[0]}_sasa_output.csv"
        SASAs = []
        TSAs = []
        volumes = []

        header = f"{os.path.splitext(fn)[0]},SASA,Total Surface Area,Volume"
        output_file.write(header)
        output_file.write("\n")
        for i, range_1 in enumerate(ranges):
            x = (rc(session, f"measure sasa {chain_format} sum {chain_format}:{range_1[0]}-{range_1[1]}"))
            SASAs.append(round(x[1], 3))
            rc(session, f"surface enclose {chain_format}:{range_1[0]}-{range_1[1]}")
            v = rc(session, f"measure volume {chain_format}:{range_1[0]}-{range_1[1]} #1.{i+1}")
            volumes.append(round(v, 3))
            t = rc(session, f"measure area {chain_format}:{range_1[0]}-{range_1[1]} #1.{i+1}")
            TSAs.append(round(t, 3))





        # Write the captured output to a text file
            output_line = f"{range_1[0]}-{range_1[1]},{str(x[1])},{str(t)},{str(v)}"

            output_file.write(output_line)
            output_file.write("\n")
        output_file.write("\n")
        output_file.write("\n")
        rc(session, "close all")


    

    # print("SASAs", SASAs)
    # print("TSAs", TSAs)
    # print("volumes", volumes)

    # Optionally, print or further process the output
    # print(f"SASA output for {fn}:\n{sasa_output}")

