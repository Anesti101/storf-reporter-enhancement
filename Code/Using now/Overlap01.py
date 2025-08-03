# -*- coding: utf-8 -*-
"""
Created on Fri Mar 21 15:42:14 2025

@author: Anesti
"""
from typing import List, Tuple, Dict
import re
import csv

# --------------------------------------------------------
# Function: parse_fasta_positions_and_sequences
# Purpose: Parse a FASTA file to extract both:
#   - Metadata from headers (like start/end positions)
#   - Full header and actual nucleotide sequence
# Logic:
#   - Group sequences by their base ID (e.g., NC_003197.2)
#   - Also track the full header and actual sequence lines
# Returns:
#   positions: {base_id: [(full_header, start, end)]}
#   sequences: {full_header: sequence}
# --------------------------------------------------------
def parse_fasta_positions_and_sequences(fasta_file: str) -> Tuple[Dict[str, List[Tuple[str, int, int]]], Dict[str, str]]:
    positions = {}
    sequences = {}
    current_header = None
    current_sequence = []

    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if current_header:
                    # Save the previous sequence before moving on
                    sequences[current_header] = ''.join(current_sequence)
                full_header = line.strip()[1:]  # Remove '>'
                match = re.search(r'^>(\S+?)_(StORF|Con-StORF)_\d+:(\d+)-(\d+)', line)
                if match:
                    seq_id = match.group(1)
                    start, end = sorted([int(match.group(3)), int(match.group(4))])
                    if seq_id not in positions:
                        positions[seq_id] = []
                    positions[seq_id].append((full_header, start, end))
                current_header = full_header
                current_sequence = []  # Reset for the next sequence
            else:
                current_sequence.append(line.strip())

        # Save last sequence
        if current_header and current_sequence:
            sequences[current_header] = ''.join(current_sequence)

    return positions, sequences

# --------------------------------------------------------
# Function: check_overlap
# Purpose: Test if two genomic ranges overlap
# Logic:
#   - Uses basic interval overlap logic
# --------------------------------------------------------
def check_overlap(range1: Tuple[int, int], range2: Tuple[int, int]) -> bool:
    return max(range1[0], range2[0]) <= min(range1[1], range2[1])

# --------------------------------------------------------
# Function: con_storf_contains_storf
# Purpose: For each Con-StORF, find overlapping StORFs
# Logic:
#   - Match by base ID
#   - Compare each con-storf range to all storfs on the same seq ID
#   - Track full headers for output
# Returns:
#   matches: {con_full_header: [storf_full_header]}
# --------------------------------------------------------
def con_storf_contains_storf(storfs: Dict[str, List[Tuple[str, int, int]]],
                              constorfs: Dict[str, List[Tuple[str, int, int]]]) -> Dict[str, List[str]]:
    contained = {}
    for seq_id in constorfs:
        if seq_id not in storfs:
            continue
        for con_full, c_start, c_end in constorfs[seq_id]:
            for storf_full, s_start, s_end in storfs[seq_id]:
                if check_overlap((c_start, c_end), (s_start, s_end)):
                    if con_full not in contained:
                        contained[con_full] = []
                    contained[con_full].append(storf_full)
    return contained

# --------------------------------------------------------
# Function: write_fasta
# Purpose: Save FASTA-formatted sequence output
# Logic:
#   - Takes a dictionary of {header: sequence}
#   - Writes them to a FASTA file
# --------------------------------------------------------
def write_fasta(output_path: str, headers: List[str], sequence_dict: Dict[str, str]):
    with open(output_path, 'w') as out:
        for header in headers:
            seq = sequence_dict.get(header)
            if seq:
                out.write(f">{header}\n")
                # Break long sequences into 60-character lines
                for i in range(0, len(seq), 60):
                    out.write(seq[i:i+60] + '\n')

# --------------------------------------------------------
# MAIN SCRIPT
# --------------------------------------------------------
if __name__ == "__main__":
    # Input files
    storf_file = "C:/Users/anest/OneDrive/Documents/CS/Y3/T2/MainP/Results/Split_500_T_Only_S/top_500_longest_storfs.fasta"
    constorf_file = "C:/Users/anest/OneDrive/Documents/CS/Y3/T2/MainP/Results/Split_500_T_Only_C/top_500_longest_Con_storfs.fasta"

    # Output overlap results
    csv_output = "C:/Users/anest/OneDrive/Documents/CS/Y3/T2/MainP/500_overlap_S_and_C.txt"

    # Output FASTA files
    storf_fasta_out = "C:/Users/anest/OneDrive/Documents/CS/Y3/T2/MainP/overlapping_storfs.fasta"
    constorf_fasta_out = "C:/Users/anest/OneDrive/Documents/CS/Y3/T2/MainP/overlapping_constorfs.fasta"

    # Parse positions and sequences from input FASTA files
    storf_coords, storf_sequences = parse_fasta_positions_and_sequences(storf_file)
    constorf_coords, constorf_sequences = parse_fasta_positions_and_sequences(constorf_file)

    # Find overlaps
    overlap_dict = con_storf_contains_storf(storf_coords, constorf_coords)

    # Write summary CSV of which Con-StORFs overlapped with which StORFs
    with open(csv_output, mode='w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Con_StORF_Full_Header", "Overlapping_StORF_Full_Headers"])
        for constorf, storfs in overlap_dict.items():
            writer.writerow([constorf, " -> ".join(storfs)])

    # Gather all headers to extract for FASTA export
    constorf_hits = list(overlap_dict.keys())
    storf_hits = [s for v in overlap_dict.values() for s in v]

    # Write overlapping sequences to FASTA files
    write_fasta(constorf_fasta_out, constorf_hits, constorf_sequences)
    write_fasta(storf_fasta_out, storf_hits, storf_sequences)
