"""
Author: Maryam KafiKang
Description: This script runs Cassiopeia analysis on phylogenetic trees to infer transmission networks.
            The process involves several key steps:
            1. Tree Processing:
               - Loads phylogenetic trees from FAVITES output
               - Extracts host labels from leaf nodes
            2. Parsimony Analysis:
               - Computes Fitch parsimony scores between all pairs of hosts
               - Generates a complete pairwise score matrix
            3. Network Inference:
               - Constructs a weighted graph from parsimony scores
               - Computes Maximum Spanning Tree (MST) to infer transmission paths
               - Preserves edge weights in the final network

Input: 
    - Newick format phylogenetic trees from FAVITES
    - Command line argument for sample index (1-500)

Output:
    - Original parsimony score matrix (output_original/*)
    - Inferred transmission network (output_infered/*)
    - Execution time measurements
"""

import os
import cassiopeia as cas
from cassiopeia.data import CassiopeiaTree
import pandas as pd
import numpy as np
import networkx as nx
import argparse
import os, time


# Set tree index
parser = argparse.ArgumentParser(description="Run Phyloscanner on FAVITES output")
parser.add_argument("i", type=int, help="Index parameter")
args = parser.parse_args()
i = args.i

# Define file paths
tree_address = f"Favites_log/FAVITES_output_logSI_contemp_T200_N100_E1_{i}/error_free_files/phylogenetic_trees/Tnet_Tree.0"
output_path_original = f"output_original/FAVITES_output_logSI_contemp_T200_N100_E1_{i}_cassiopeia.txt"
output_path_infered = f"output_infered/FAVITES_output_logSI_contemp_T200_N100_E1_{i}_cassiopeia.txt"

# Ensure output directories exist
os.makedirs(os.path.dirname(output_path_original), exist_ok=True)
os.makedirs(os.path.dirname(output_path_infered), exist_ok=True)

# Check if the tree file exists before proceeding
if not os.path.exists(tree_address):
    print(f"Error: Tree file '{tree_address}' not found.")
    exit()

# Start timer
start_time = time.time()

# Load the phylogenetic tree
tree = CassiopeiaTree()
tree.populate_tree(tree_address)

# Compute the Fitch Parsimony Score
leaves = tree.leaves
host_labels = [leaf.split("_")[0] for leaf in leaves]
tree.cell_meta = pd.DataFrame(index=leaves, data={"hosts": host_labels})

# Compute the parsimony score matrix
parsimony_score_df = cas.tl.fitch_count(tree, meta_item="hosts")

# End timer
execution_time = time.time() - start_time

# Save Parsimony Score Matrix
with open(output_path_original, 'w') as f:
    f.write(parsimony_score_df.to_string(index=True))  # Convert DataFrame to string
    f.write(f"\nRunning Time (seconds): {execution_time:.5f}\n")

# Construct graph edges from parsimony score matrix
edges = [
    (i, j, parsimony_score_df.loc[i, j])
    for i in parsimony_score_df.index
    for j in parsimony_score_df.columns
    if i != j  # Exclude self-loops
]

# Create an undirected graph
G = nx.Graph()
G.add_weighted_edges_from(edges)

# Compute the Maximum Spanning Tree (MST)
mst = nx.maximum_spanning_tree(G, weight="weight")

# Convert MST to DataFrame
mst_edges = pd.DataFrame(mst.edges(data=True), columns=["Node1", "Node2", "Attributes"])
mst_edges["Weight"] = mst_edges["Attributes"].apply(lambda x: x["weight"])
mst_edges.drop(columns=["Attributes"], inplace=True)

# Save MST edges to file
mst_edges.to_csv(output_path_infered, sep="\t", index=False, header=False)

# Append execution time to the output file
with open(output_path_infered, 'a') as f:
    f.write(f"\nRunning Time (seconds): {execution_time:.5f}\n")

# Print the final results
print("\nFitch Parsimony Score Matrix:")
print(parsimony_score_df)
print("\nMaximum Spanning Tree (MST) Edges:")
print(mst_edges)
print(f"\nProcessing completed in {execution_time:.5f} seconds.")


        
