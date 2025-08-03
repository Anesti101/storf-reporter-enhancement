# -*- coding: utf-8 -*-
"""
Created on Wed Apr  2 16:20:04 2025

@author: Anesti
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# --------------------------------------------------------
# Function: process_blast_files
# Purpose: Load two BLAST tabular files and annotate their metrics
# Logic:
#   - Load both StORF and Con-StORF BLAST result files into DataFrames
#   - Label them with their type and merge into one combined DataFrame
#   - Calculate query and subject length, and coverage percentages
#   - Classify each row by percent identity, E-value, bit score, and coverage
#   - Return the annotated DataFrame for downstream plotting
# --------------------------------------------------------
def process_blast_files(storf_file_path, constorf_file_path):
    # Define columns for BLAST format 6
    columns = [
        "Query_ID", "Subject_ID", "Percent_Identity", "Alignment_Length",
        "Mismatches", "Gap_Openings", "Q_Start", "Q_End",
        "S_Start", "S_End", "E_Value", "Bit_Score"
    ]

    # Load the files
    storf_df = pd.read_csv(storf_file_path, sep='\t', header=None, names=columns)
    constorf_df = pd.read_csv(constorf_file_path, sep='\t', header=None, names=columns)

    # Label them
    storf_df["Type"] = "StORF"
    constorf_df["Type"] = "Con-StORF"

    # Combine
    df = pd.concat([storf_df, constorf_df], ignore_index=True)

    # Compute additional fields
    df["Query_Length"] = abs(df["Q_End"] - df["Q_Start"]) + 1
    df["Subject_Length"] = abs(df["S_End"] - df["S_Start"]) + 1
    df["Query_Coverage_%"] = df["Alignment_Length"] / df["Query_Length"] * 100
    df["Subject_Coverage_%"] = df["Alignment_Length"] / df["Subject_Length"] * 100

    # ------------------ CATEGORISATION FUNCTIONS ------------------

    # Categorise percent identity into descriptive strength bins
    def identity_category(val):
        if val >= 90:
            return "Strong (≥90%)"
        elif 80 <= val < 90:
            return "Moderate (80–90%)"
        elif val < 50:
            return "Weak (<50%)"
        return "Other"

    df["Identity_Strength"] = df["Percent_Identity"].apply(identity_category)

    # Categorise E-values into meaningful significance levels
    def evalue_category(val):
        try:
            val = float(val)
            if val <= 1e-10:
                return "Very Strong (≤1e-10)"
            elif val <= 1e-5:
                return "Strong (≤1e-5)"
            elif val > 1e-3:
                return "Weak (>1e-3)"
        except:
            return "Unknown"
        return "Other"

    df["EValue_Strength"] = df["E_Value"].apply(evalue_category)

    # Categorise bit scores to reflect match quality
    def bitscore_category(val):
        if val > 200:
            return "Very Strong (>200)"
        elif val > 100:
            return "Strong (>100)"
        elif val < 50:
            return "Weak (<50)"
        return "Other"

    df["BitScore_Strength"] = df["Bit_Score"].apply(bitscore_category)

    # Categorise matches based on both query and subject coverage
    def coverage_category(row):
        q_cov = row["Query_Coverage_%"]
        s_cov = row["Subject_Coverage_%"]
        if q_cov > 80 and s_cov > 50:
            return "Strong Coverage"
        elif q_cov < 30 or s_cov < 30:
            return "Weak Coverage"
        return "Moderate Coverage"

    df["Coverage"] = df.apply(coverage_category, axis=1)
    
    return df


# --------------------------------------------------------
# Function: plot_comparison
# Purpose: Plot and compare the distribution of a categorical metric by sequence type
# Logic:
#   - Uses seaborn countplot to show frequency of each category
#   - Groups counts by the selected category and Type (StORF or Con-StORF)
#   - Returns a DataFrame with the grouped counts
# --------------------------------------------------------
def plot_comparison(df, category_column, title):
    plt.figure(figsize=(10, 6))
    sns.countplot(data=df, x=category_column, hue="Type", palette="Set2")
    plt.title(title)
    plt.xlabel("Category")
    plt.ylabel("Number of Hits")
    plt.xticks(rotation=15)
    plt.legend(title="Sequence Type")
    plt.tight_layout()
    plt.show()
    return df[[category_column, "Type"]].value_counts().reset_index(name="Count")


# --------------------------------------------------------
# MAIN EXECUTION: Run the workflow on specific input files
# Logic:
#   - Load and annotate both StORF and Con-StORF BLAST hits
#   - Visualise the distribution of key metrics via bar plots
#   - Store the count results for further analysis/reporting
# --------------------------------------------------------
storf_file_path = "C:/Users/anest/OneDrive/Documents/CS/Y3/T2/MainP/Results/Split_500_T_Only_S/YUKKS3SW013-Alignment-HitTable.csv" 
constorf_file_path = "C:/Users/anest/OneDrive/Documents/CS/Y3/T2/MainP/Results/Split_500_T_Only_C/YUM5ZGK7016-Alignment-HitTable.csv" 
df_combined = process_blast_files(storf_file_path, constorf_file_path)

# Show all plot categories
identity_counts = plot_comparison(df_combined, "Identity_Strength", "Percent Identity Categories by Type")
evalue_counts = plot_comparison(df_combined, "EValue_Strength", "E-Value Categories by Type")
bitscore_counts = plot_comparison(df_combined, "BitScore_Strength", "Bit Score Categories by Type")
coverage_counts = plot_comparison(df_combined, "Coverage", "Coverage Categories by Type")
