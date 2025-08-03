# -*- coding: utf-8 -*-
"""
Created on Sun Feb 23 16:21:16 2025

@author: Anesti
"""

import os 


""""Relabels the sequences form the files and merging them"""
def relabel_and_merge(storf_file, con_storf_file, output_file):
    """It reads the files and relabels the sequences form a Fasta file with a the given prefix """
    def process_fasta(input_file, prefix):
        sequences = []
        with open(input_file, 'r') as infile:
        # initailize counter for seq numbering 
            counter = 1 
            for line in infile:
                if line.startswith(">"): 
                    # Replaces the header with the prefix and counter
                    sequences.append(f">{prefix}{counter}\n")
                    # Increment conter for the next seq
                    counter += 1
                else:
                    # for seq line add them as they are
                    sequences.append(line)
        return sequences
    
    # check if both input files exist before processing them
    if not os.path.exists(storf_file) or not os.path.exists(con_storf_file):
        print("error: one or both input files not found.")
        return
    
    # gives the prefix S to the Storf file and prefix C for the Con-Storf
    storf_sequences = process_fasta(storf_file, "S")
    con_storf_sequences = process_fasta(con_storf_file, "C")
    
    # Merge them and write them to the output file
    with open(output_file, "w") as outfile:
        outfile.writelines(storf_sequences + con_storf_sequences)
    print(f"Merged file '{output_file}' created succesfully")

# File paths 
storf_file = "C:/Users/anest/OneDrive/Documents/CS/Y3/T2/MainP/Results/Only_S/S_GCF_000006945_StORF-Finder.fasta"  
con_storf_file = "C:/Users/anest/OneDrive/Documents/CS/Y3/T2/MainP/Results/Only_C/C_GCF_000006945_StORF-Finder.fasta"  
output_file = "C:/Users/anest/OneDrive/Documents/CS/Y3/T2/MainP/Results/merged_storf_con_storf.fasta"  

# Run the merge function 
relabel_and_merge(storf_file, con_storf_file, output_file)
