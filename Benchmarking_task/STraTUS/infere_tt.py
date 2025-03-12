"""
Author: Maryam KafiKang
Description: A companion script that processes STraTUS output (_tt.csv files) to infer 
                    transmission trees. It works by:
                    1. Extracting edges and their weights from STraTUS CSV output
                    2. Building a weighted graph from these edges
                    3. Computing a maximum spanning tree to determine the most likely
                       transmission paths
                    4. Saving the results in a format suitable for evaluation
"""

import pandas as pd
import networkx as nx
import os
from collections import Counter

def extract_edges(csv_file):
    """Extract edges and compute weights based on repetition."""
    df = pd.read_csv(csv_file)
    edge_list = []
    
    # Iterate through each column (excluding the 'child' column)
    for col in df.columns[1:]:
        for child, parent in zip(df['child'], df[col]):
            if pd.notna(parent) and parent != "root":  # Ignore 'root' entries
                edge_list.append((parent, child))
    
    # Count occurrences of each edge
    edge_weights = Counter(edge_list)
    return edge_weights

def compute_max_spanning_tree(edge_weights, output_file):
    """Constructs a weighted graph, computes the Maximum Spanning Tree, and saves it."""
    G = nx.Graph()
    
    # Add weighted edges to graph
    for (parent, child), weight in edge_weights.items():
        G.add_edge(parent, child, weight=weight)
    
    # Compute Maximum Spanning Tree
    mst = nx.maximum_spanning_tree(G, weight='weight')
    
    # Save MST to file
    with open(output_file, "w") as f:
        for u, v, data in mst.edges(data=True):
            weight = data.get('weight', 0)  # Extract the weight, default to 0 if missing
            f.write(f"{u.split('H')[0]} {v.split('H')[0]}   {weight} \n")
    
    print(f"Maximum spanning tree saved to {output_file}")

if __name__ == "__main__":
    for i in range(1, 501):
        output_dir = 'Inferred_tt'
        input_csv = f"stratus_output/FAVITES_output_logSI_contemp_T200_N100_E1_{i}/_tt.csv"
        if os.path.exists(input_csv):
            os.makedirs(output_dir, exist_ok=True)
            output_txt = f"{output_dir}/FAVITES_output_logSI_contemp_T200_N100_E1_{i}_tt.txt"
            #if not os.path.exists(output_txt):
            edges = extract_edges(input_csv)
            compute_max_spanning_tree(edges, output_txt)
            
