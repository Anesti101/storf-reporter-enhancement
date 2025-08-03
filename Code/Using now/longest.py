"""
Created on Sun Feb 16 16:43:19 2025

@author: Anesti
"""

import pandas as pd

# Set path to the raw BLAST results (tab-separated values)
blast_results_file = 'C:/Users/anest/OneDrive/Documents/CS/Y3/T2/MainP/Results/S_and_C/S_and_C_blastn_results.tsv'

# Set path to save the output Excel file after comparison
output_excel_file = 'C:/Users/anest/OneDrive/Documents/CS/Y3/T2/MainP/Results/S_and_C_length_analysis.xlsx'

# Define the expected BLAST column headers
columns = [
    'Query_ID', 'Subject_ID', 'Percent_Identity', 'Alignment_Length',
    'Mismatches', 'Gap_Openings', 'Q_Start', 'Q_End', 'S_Start',
    'S_End', 'E_Value', 'Bit_Score'
]

# Read the BLAST results into a pandas DataFrame
blast_results = pd.read_csv(blast_results_file, sep='\t', header=None, names=columns)


# For each Subject_ID find the row with the longest alignment length
longest_storf = blast_results.loc[blast_results.groupby('Subject_ID')['Alignment_Length'].idxmax()]


# For each Subject_ID, find the row with the lowest E-value (best BLAST hit)
best_match_storf = blast_results.loc[blast_results.groupby('Subject_ID')['E_Value'].idxmin()]


# Merge both longest and best match DataFrames on Subject_ID
# Add suffixes to distinguish between the two sources
comparison_df = longest_storf.merge(best_match_storf, on='Subject_ID', suffixes=('_Longest', '_Best'))


# Create a boolean column to check whether the longest STORF is also the best match
comparison_df['Is_Longest_Best'] = comparison_df['Query_ID_Longest'] == comparison_df['Query_ID_Best']


# Write the comparison results to an Excel file for further analysis
with pd.ExcelWriter(output_excel_file, engine='openpyxl') as writer:
    comparison_df.to_excel(writer, sheet_name='Longest_vs_Best', index=False)

# Calculate how often the longest STORF is also the best match
longest_selected = comparison_df['Is_Longest_Best'].mean() * 100
print(f"The longest STORF is the best match in selected: {longest_selected:.2f}% of cases.")

# Display key comparison columns for inspection
print(comparison_df[['Query_ID_Longest', 'Query_ID_Best', 'Subject_ID',
                     'Percent_Identity_Longest', 'E_Value_Longest',
                     'E_Value_Best', 'Is_Longest_Best']])