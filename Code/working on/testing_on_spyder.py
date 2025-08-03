# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 13:04:16 2025

@author: Anesti
"""

import os
#import sys
#sys.path.append(r"C:\Users\anest\AppData\Roaming\Python\Python38\site-packages")
#import StORF_Reporter
#print("StORF-Reporter is now accessible!")

# Define file paths
output_dir = "C:/Users/anest/Downloads/Ang76_Part_1/results"
fasta_file = "C:/Users/anest/Downloads/ncbi_dataset/ncbi_dataset/data/GCF_000006945.2/GCF_000006945.2_ASM694v2_genomic.fna"
gff_file = "C:/Users/anest/Downloads/ncbi_dataset/ncbi_dataset/data/GCF_000006945.2/genomic.gff"

# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)
print("Output directory is ready:", output_dir)

# Step 1: Extract Unannotated Regions (URs)
print("Running UR-Extractor...")
os.system(f'UR-Extractor -f "{fasta_file}" -gff "{gff_file}" -odir "{output_dir}" -gz True')

# Check if UR files are created
ur_fasta = os.path.join(output_dir, "GCF_000006945_UR.fasta.gz")
if os.path.exists(ur_fasta):
    print(f"UR file created: {ur_fasta}")
else:
    print("UR file not created. Please check the UR-Extractor step.")

# Step 2: Find StORFs in Unannotated Regions
if os.path.exists(ur_fasta):
    print("Running StORF-Finder...")
    os.system(f'StORF-Finder -f "{ur_fasta}" -odir "{output_dir}/GCF_000006945_StORFs" -aa True -gff True -gz True')
else:
    print("Skipping StORF-Finder since UR file is missing.")

# Step 3: Run StORF-Reporter to Automate the Full Process
print("Running StORF-Reporter...")
os.system(f'StORF-Reporter -anno Pyrodigal -p "{fasta_file}" -odir "{output_dir}" -aa True -sout True')
