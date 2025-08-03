# -*- coding: utf-8 -*-
"""
Created on Sat Mar 29 15:29:51 2025

@author: Anesti
"""
from Bio import SeqIO
import os
# --------------------------------------------------------
# Function: split_longest_sequences
# Purpose: Split a FASTA file into:
#   - Top N longest sequences
#   - The remaining sequences
# Logic:
#   - Load all sequences from input FASTA
#   - Sort sequences by length (ascending)
#   - Take top N longest and the rest
#   - Write each group to separate FASTA files
# Input:
#   - input_fasta: path to FASTA file
#   - output_dir: folder where outputs will be saved
#   - top_n: number of longest sequences to extract
# Output:
#   - Two FASTA files: top_N_longest_storfs.fasta, remaining_storfs.fasta
# --------------------------------------------------------
def split_longest_sequences(input_fasta, output_dir, top_n=20):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Load all sequences from the FASTA file
    records = list(SeqIO.parse(input_fasta, "fasta"))

    # Sort sequences by length ascending 
    sorted_records = sorted(records, key=lambda r: len(r.seq), reverse=False) #

    # Split into top N longest and the rest
    longest_seqs = sorted_records[:top_n] 
    other_seqs = sorted_records[top_n:]

    # Output file paths
    longest_path = os.path.join(output_dir, f"top_{top_n}_longest_storfs.fasta") 
    others_path = os.path.join(output_dir, f"remaining_storfs.fasta") # wri

    # wrtie the sequences to the output files
    SeqIO.write(longest_seqs, longest_path, "fasta") # write the longest sequences 
    SeqIO.write(other_seqs, others_path, "fasta") # write the remaining sequences

    # print cofiramation messages
    print(f"Saved top {top_n} longest sequences to: {longest_path}")
    print(f"Saved remaining {len(other_seqs)} sequences to: {others_path}")

# Usage 
input_fasta_path = "C:/Users/anest/OneDrive/Documents/CS/Y3/T2/MainP/Results/S_and_C/ResultsGCF_000006945_StORF-Finder.fasta"
output_directory = "C:/Users/anest/OneDrive/Documents/CS/Y3/T2/MainP/Results/Split_500_B"
top_n = 500  # Number of longest sequences to extract

split_longest_sequences(input_fasta_path, output_directory, top_n)
