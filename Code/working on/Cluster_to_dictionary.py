# -*- coding: utf-8 -*-
"""
Created on Sun Feb 23 16:54:06 2025

@author: Anesti
"""
import os 
import re

""" Parse the clustered data from a file. Extracts cluster names and sequence details."""
def parse_clustered_file(file_path):
    
    # Initialize an empty dictionary to store clusters and their associated sequence details.
    clusters = {}
    # This variable will hold the identifier of the current cluster being processed.
    current_cluster = None 
    
    # Check if the file exists
    if not os.path.exists(file_path): 
        print(f"Error: File '{file_path}' not found.")
        return clusters
    
    # Read file line by line
    with open(file_path, 'r') as infile: 
        for line in infile:
            # Remove any leading and trailing whitespace from the line.
            line = line.strip() 
            # If the line indicates the start of a new cluster
            if line.startswith(">Cluster"):
                # Split the line and use the second token as the cluster identifier
                current_cluster = line.split()[1] # New cluster detected 
                # Initialize an empty list in the dictionary for this new cluster.
                clusters[current_cluster] = []
            else:
                # For lines that do not start a new cluster, extract seq details using a regular expression
                match = re.match(r"(\d+)\s+(\d+)nt,\s+>(S|C)\d+\.\.\.\s+at\s+[-+]/(\d+\.\d+)%?", line)
                if match:
                    seq_id = match.group(1)
                    length = int(match.group(2))
                    origin = match.group(3)
                    similarity = float(match.group(4))
                    # Append a tuple containing the extracted details to the list for the current cluster
                    clusters[current_cluster].append((seq_id, length, origin, similarity))
    return clusters
