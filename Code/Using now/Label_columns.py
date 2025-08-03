# -*- coding: utf-8 -*-
"""
Created on Sat Mar 29 18:17:55 2025

@author: Anesti
"""
# --------------------------------------------------------
# Script: Label BLAST Hit Table
# Purpose: Clean and label a BLAST alignment hit table
# Logic:
#   - Load raw BLAST CSV file with no column headers
#   - Assign proper column names based on standard BLAST format
#   - Save the cleaned and labelled data to a new CSV file
# Input:
#   - input_file: raw BLAST output (CSV format, comma-separated)
# Output:
#   - output_file: labelled version of the same BLAST table
# --------------------------------------------------------
import pandas as pd

input_file = "C:/Users/anest/OneDrive/Documents/CS/Y3/T2/MainP/Results/Top_n_BLAST/Shortest/B_500_YRV7YV5N016-Alignment-HitTable.csv"  
output_file = "C:/Users/anest/OneDrive/Documents/CS/Y3/T2/MainP/Results/Top_n_BLAST/Shortest/B_500_L_results_labeled.csv"

# Read using comma as separator
df = pd.read_csv(input_file, sep=",", header=None)

print(df.head())
print(df.shape)

# Assign BLAST column headers
columns = [
    "Query_ID", "Subject_ID", "Percent_Identity", "Alignment_Length",
    "Mismatches", "Gap_Openings", "Q_Start", "Q_End",
    "S_Start", "S_End", "E_Value", "Bit_Score"
]

df.columns = columns

# Save cleaned and labelled file
df.to_csv(output_file, index=False)
print(f"Labeled file saved as: {output_file}")