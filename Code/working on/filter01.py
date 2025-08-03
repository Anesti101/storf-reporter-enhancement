# -*- coding: utf-8 -*-
"""
Created on Sat Mar 22 17:37:55 2025

@author: Anesti
"""

import pandas as pd
import re
from typing import List, Tuple, Dict

def parse_fasta_positions(fasta_file: str) -> Dict[str, Tuple[int, int]]:
    """Extract coordinates from FASTA headers into a dict of ID: (start, end)."""
    coords = {}
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                match = re.search(r'^>(\S+?)_(StORF|Con-StORF)_\d+:(\d+)-(\d+)', line)
                if match:
                    id_full = match.group(0)[1:].split(';')[0].strip()  # clean full ID
                    start, end = sorted([int(match.group(3)), int(match.group(4))])
                    coords[id_full] = (start, end)
    return coords

def check_overlap(range1: Tuple[int, int], range2: Tuple[int, int]) -> bool:
    return max(range1[0], range2[0]) <= min(range1[1], range2[1])

def calculate_overlap_length(range1: Tuple[int, int], range2: Tuple[int, int]) -> int:
    """Calculate the length of the overlapping region between two ranges."""
    return max(0, min(range1[1], range2[1]) - max(range1[0], range2[0]) + 1)

def analyse_top_overlap(blast_file: str, storf_file: str, constorf_file: str, output_file: str, top_n: int = 100):
    # Load coordinates
    storf_coords = parse_fasta_positions(storf_file)
    constorf_coords = parse_fasta_positions(constorf_file)

    # Load blast result (format 6 assumed)
    cols = ["Query_ID", "Subject_ID", "%_Identity", "Alignment_Length", "Mismatch", "Gap_Openings",
            "Q_Start", "Q_End", "S_Start", "S_End", "E_Value", "Bit_Score"]
    blast_df = pd.read_csv(blast_file, sep='\t', names=cols)

    # Clean and simplify Query_IDs
    blast_df['Clean_Query_ID'] = blast_df['Query_ID'].str.strip().str.split(';').str[0]

    # Get top N unique query IDs
    top_queries_df = blast_df.drop_duplicates(subset='Clean_Query_ID').head(top_n)
    top_queries = top_queries_df['Clean_Query_ID'].tolist()

    print("Parsed Con-StORF IDs:", list(constorf_coords.keys())[:5])
    print("Top BLAST Query IDs:", top_queries[:5])

    con_storf_hits = [qid for qid in top_queries if 'Con-StORF' in qid]
    results = []
    overlap_count = 0
    storf_selected = 0

    for con_id in con_storf_hits:
        con_id_clean = con_id.strip()
        con_range = constorf_coords.get(con_id_clean)

        # If direct match doesn't work, try partial match
        if not con_range:
            for key in constorf_coords:
                if con_id_clean in key:
                    con_range = constorf_coords[key]
                    con_id_clean = key
                    break

        if not con_range:
            continue

        # Search for overlapping StORFs
        for storf_id, storf_range in storf_coords.items():
            if check_overlap(con_range, storf_range):
                overlap_count += 1
                storf_rank = top_queries.index(storf_id) if storf_id in top_queries else -1
                con_rank = top_queries.index(con_id)
                blast_row = blast_df[blast_df['Clean_Query_ID'] == con_id_clean].iloc[0]

                pick = "Con-StORF"
                if storf_rank != -1 and (con_rank - storf_rank) >= 20:
                    pick = "StORF"
                    storf_selected += 1

                overlap_len = calculate_overlap_length(con_range, storf_range)
                storf_len = storf_range[1] - storf_range[0] + 1
                constorf_len = con_range[1] - con_range[0] + 1

                results.append({
                    'Con-StORF': con_id_clean,
                    'StORF': storf_id,
                    'Con_StORF_Rank': con_rank,
                    'StORF_Rank': storf_rank,
                    'StORF_Higher': storf_rank != -1 and storf_rank < con_rank,
                    'Preferred_Pick': pick,
                    'StORF_Length': storf_len,
                    'Con_StORF_Length': constorf_len,
                    'Overlap_Length': overlap_len,
                    '%_Identity': blast_row['%_Identity'],
                    'Bit_Score': blast_row['Bit_Score'],
                    'E_Value': blast_row['E_Value'],
                    'Alignment_Length': blast_row['Alignment_Length']
                })
                break

    # Save overlap summary to a text file
    with open("ecoli_overlap_summary.txt", "w") as f:
        f.write(f"Total overlaps detected: {overlap_count}\n")
        f.write(f"StORF selected over Con-StORF: {storf_selected} times\n")
        f.write(f"Con-StORF selected: {overlap_count - storf_selected} times\n")

    print(f"Summary written to ecoli_overlap_summary.txt")

    results_df = pd.DataFrame(results)
    results_df.to_csv(output_file, index=False)
    print(f"Results written to {output_file}")
    return results_df

if __name__ == "__main__":
    
    storf_file = "C:/Users/anest/OneDrive/Documents/CS/Y3/T2/MainP/Results/Only_S/S_GCF_000006945_StORF-Finder.fasta"         
    constorf_file = "C:/Users/anest/OneDrive/Documents/CS/Y3/T2/MainP/Results/Only_C/C_GCF_000006945_StORF-Finder.fasta"  
    output_file = "C:/Users/anest/OneDrive/Documents/CS/Y3/T2/MainP/Results/Answer8 .csv"     
    blast_file = "C:/Users/anest/OneDrive/Documents/CS/Y3/T2/MainP/s_and_c_Blast1.txt"

    analyse_top_overlap(blast_file, storf_file, constorf_file, output_file, top_n=500)
