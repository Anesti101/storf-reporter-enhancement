# -*- coding: utf-8 -*-
"""
Created on Sat Mar 29 18:45:07 2025

@author: Anesti
"""
import pandas as pd
import matplotlib.pyplot as plt

# --------------------------------------------------------
# Function: load_blast_file
# Purpose: Load a BLAST CSV file into a DataFrame with metrics
# Logic:
#   - Open the file line by line and split by commas
#   - Validate that each line has 12 columns (standard BLAST output)
#   - Convert numeric fields to appropriate types (int or float)
#   - Append data to a dictionary
#   - Add a column indicating whether it's StORF or Con-StORF
#   - Calculate query length and query coverage
# Input:
#   - file_path: path to the CSV file
#   - sequence_type: string label (e.g. 'StORF' or 'Con-StORF')
# Output:
#   - DataFrame containing cleaned and labelled data
# --------------------------------------------------------
def load_blast_file(file_path, sequence_type):


    print(f" Loading: {file_path} as {sequence_type}")
    
    # Define dictionary to hold columns
    data = {
        "Query_ID": [], "Subject_ID": [], "Percent_Identity": [], "Alignment_Length": [],
        "Mismatches": [], "Gap_Openings": [], "Q_Start": [], "Q_End": [],
        "S_Start": [], "S_End": [], "E_Value": [], "Bit_Score": [], "Type": []
    }

    # Read the file line by line and split into columns
    with open(file_path, 'r') as f:
        for i, line in enumerate(f, start=1):
            parts = line.strip().split(',')
            if len(parts) != 12:
                print(f" Skipping line {i}: not 12 columns")
                continue

            try:
                # Append parsed values to corresponding columns
                data["Query_ID"].append(parts[0])
                data["Subject_ID"].append(parts[1])
                data["Percent_Identity"].append(float(parts[2]))
                data["Alignment_Length"].append(int(parts[3]))
                data["Mismatches"].append(int(parts[4]))
                data["Gap_Openings"].append(int(parts[5]))
                data["Q_Start"].append(int(parts[6]))
                data["Q_End"].append(int(parts[7]))
                data["S_Start"].append(int(parts[8]))
                data["S_End"].append(int(parts[9]))
                data["E_Value"].append(float(parts[10]) if parts[10].lower() != "nan" else None)
                data["Bit_Score"].append(float(parts[11]))
                data["Type"].append(sequence_type)
            except ValueError as e:
                print(f" Error parsing line {i}: {e}")
                continue

    # Convert dictionary to DataFrame
    df = pd.DataFrame(data)
    print(f" {sequence_type} loaded: {len(df)} valid rows")

    # Compute query length and coverage
    df["Query_Length"] = abs(df["Q_End"] - df["Q_Start"]) + 1
    df["Query_Coverage"] = df["Alignment_Length"] / df["Query_Length"] * 100

    return df


# --------------------------------------------------------
# Function: plot_metrics
# Purpose: Generate and save boxplots for BLAST metrics
# Logic:
#   - Take a DataFrame with BLAST results
#   - Plot 3 boxplots:
#       • Percent Identity
#       • Bit Score
#       • Query Coverage
#   - Save the figure as a PNG image and also display it
# Input:
#   - df: BLAST result DataFrame
#   - sequence_type: label used in figure title and filename
# Output:
#   - Boxplot image saved to disk
# --------------------------------------------------------
def plot_metrics(df, sequence_type):

    metrics = ["Percent_Identity", "Bit_Score", "Query_Coverage"]

    # Create subplots
    fig, axs = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle(f"Metrics Distribution for {sequence_type}", fontsize=16)

    # Plot each metric in a boxplot
    for i, metric in enumerate(metrics):
        axs[i].boxplot(df[metric].dropna(), vert=True)
        axs[i].set_title(metric)
        axs[i].set_ylabel("Value")
        axs[i].grid(True, linestyle='--', alpha=0.5)

    # Save and display the figure
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(f"{sequence_type}_metrics_boxplots.png", dpi=300)
    plt.show()



# File paths for your BLAST results 
storf_file = "C:/Users/anest/OneDrive/Documents/CS/Y3/T2/MainP/S_YUWCUTNY016-Alignment-HitTable.csv"
constorf_file = "C:/Users/anest/OneDrive/Documents/CS/Y3/T2/MainP/Con_YUWMX1WN013-Alignment-HitTable.csv"

# Run pipeline for both types 
# Load and parse BLAST results for StORFs
storf_df = load_blast_file(storf_file, "StORF")

# Load and parse BLAST results for Con-StORFs
constorf_df = load_blast_file(constorf_file, "Con-StORF")

# --------------------------------------------------------
# IF storf_df is not empty:
#    CALL plot_metrics on storf_df
#
# IF constorf_df is not empty:
#    CALL plot_metrics on constorf_df#
#
# ELSE:
#   Print message indicating no data to plot
# --------------------------------------------------------

if not storf_df.empty:                   
    plot_metrics(storf_df, "StORF")
else:
    print(" Skipping plot for StORF — empty DataFrame.")

if not constorf_df.empty:
    plot_metrics(constorf_df, "Con-StORF")
else:
    print(" Skipping plot for Con-StORF — empty DataFrame.")



