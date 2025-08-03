# -*- coding: utf-8 -*-
"""
Created on Tue Apr  1 15:52:32 2025

@author: Anesti
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as mticker  # for clean y-axis tick formatting


longest_file = "C:/Users/anest/OneDrive/Documents/CS/Y3/T2/MainP/Results/Top_n_BLAST/Longest/100_YJPS09WP013-Alignment-HitTable.csv"

# Define column names based on BLAST output format 6
columns = [
    "Query_ID", "Subject_ID", "Percent_Identity", "Alignment_Length",
    "Mismatches", "Gap_Openings", "Q_Start", "Q_End",
    "S_Start", "S_End", "E_Value", "Bit_Score"
]

# Read longest results into a DataFrame and label
long_df = pd.read_csv(longest_file, header=None, names=columns, dtype=str)
long_df["Length_Group"] = "Longest"

# Convert numeric columns after loading
numeric_cols = ["Percent_Identity", "Alignment_Length", "Mismatches", "Gap_Openings",
                "Q_Start", "Q_End", "S_Start", "S_End", "E_Value", "Bit_Score"]

for col in numeric_cols:
    long_df[col] = pd.to_numeric(long_df[col], errors="coerce")


# === Classify StORF Type based on metadata in Query_ID ===

# Define function to extract 'StORF_Type' from the metadata string in Query_ID
def classify_by_metadata(qid):
    qid = str(qid)  # Ensure it's a string
    if "StORF_Type=" in qid:
        try:
            # Split at 'StORF_Type=' and take the value before the next ';'
            return qid.split("StORF_Type=")[1].split(";")[0]
        except IndexError:
            return "Unknown"  # Catch malformed entries
    return "Unknown"  # If key is missing

# Apply the classification function to the Query_ID column
long_df["Type"] = long_df["Query_ID"].apply(classify_by_metadata)
long_df["Type"] = long_df["Type"].astype(str)  # Ensure type is string

# Set Seaborn plot style
sns.set(style="whitegrid", context="talk")

# === Plot 1: Percent Identity comparison ===
plt.figure(figsize=(12, 6)) # set figure size
sns.boxplot(
    data=long_df,
    x="Length_Group",
    y="Percent_Identity",
    hue="Type",  # Colour by StORF/con-StORF
    palette="Set2"
)


# Add labels and title
plt.title("Percent Identity of Shortest vs Longest Sequences", fontsize=16)
plt.ylabel("Percent Identity (%)", fontsize=12)
plt.xlabel("Length Group", fontsize=12)

# Format y-axis ticks: step every 5%, and show as "XX%"
plt.gca().yaxis.set_major_locator(mticker.MultipleLocator(10))
plt.gca().yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f'{x:.0f}%'))

# Set font size 
plt.yticks(fontsize=10)
plt.xticks(fontsize=12)

# Add grid for clarity
plt.grid(axis='y', linestyle='--', alpha=0.3)
plt.legend(title="Sequence Type") # Add legend
plt.tight_layout() # Ensure layout fits
plt.savefig("longest_percent_identity_boxplot.png", dpi=300)
plt.show()



# === Plot 2: Count of StORF vs con-StORF ===
plt.figure(figsize=(8, 6)) # Set figure size
sns.countplot(
    data=long_df,
    x="Length_Group",
    hue="Type",  # Count bars by StORF type
    palette="Set2"
)

# Add labels and title
plt.title("Count of StORF vs con-StORF in Longest Sequences", fontsize=16)
plt.xlabel("Length Group", fontsize=12)
plt.ylabel("Number of Hits", fontsize=12)

# Set font sizes
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.grid(axis='y', linestyle='--', alpha=0.3) # Add light grid
plt.legend(title="Sequence Type") 
plt.tight_layout()
plt.savefig("longest_storf_type_count.png", dpi=300)
plt.show()


