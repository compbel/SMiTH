"""
Author: Maryam KafiKang
Description: This script infers transmission trees from TNet network output
            by computing maximum spanning trees from weighted edges.
            It processes multiple simulation outputs and preserves execution times.
"""

import networkx as nx
import time
import os

def read_graph_from_file(filename):
    """
    Reads a text file containing weighted edges and extracts them into a graph.
    """
    G = nx.Graph()  # Create an undirected graph
    execution_time = None  # Store execution time if available
    with open(filename, "r") as file:
        for line in file:
            line = line.strip()
            if not line:
                continue  # Skip empty lines
            
            # Check for execution time in the last line
            if "Running Time" in line:
                try:
                    execution_time = float(line.split(":")[-1].strip())  # Extract time
                except ValueError:
                    execution_time = None  # Handle any parsing errors
                continue  # Skip adding this line to the graph

            # Extract edge and weight
            parts = line.split()
            if len(parts) != 2:
                continue  # Skip malformed lines
            
            edge, weight = parts[0], parts[1]
            weight = int(weight)  # Convert weight to integer
            
            # Extract nodes (assuming format "A->B")
            if "->" in edge:
                node1, node2 = edge.split("->")
                G.add_edge(node1, node2, weight=weight)  # Add edge to graph

    return G, execution_time

def compute_maximum_spanning_tree(G):
    """
    Computes the Maximum Spanning Tree (MST) of the given graph.
    """
    return nx.maximum_spanning_tree(G, weight='weight')

def write_mst_to_file(MST, output_filename, execution_time):
    """
    Writes the edges of the Maximum Spanning Tree to an output file and includes execution time.
    """
    with open(output_filename, "w") as file:
        for u, v, data in MST.edges(data=True):
            file.write(f"{u} {v} \n")
    
        if execution_time is not None:
            file.write(f"\nRunning Time (seconds): {execution_time}\n")

    print(f"Maximum Spanning Tree saved to {output_filename}")

# Example usage
for i in range(1,501):
    input_file = f"tnet_output_updated/FAVITES_output_logSI_contemp_T200_N100_E1_{i}_TNet_Network.txt"  # Replace with your input file
    output_file = f"inferedTree/FAVITES_output_logSI_contemp_T200_N100_E1_{i}_TNet_Network.txt"

    # Read graph from file
    if os.path.exists(input_file):
        graph, recorded_execution_time = read_graph_from_file(input_file)

        # Compute Maximum Spanning Tree
        mst = compute_maximum_spanning_tree(graph)


        # Write MST to file with execution time
        write_mst_to_file(mst, output_file, recorded_execution_time)
