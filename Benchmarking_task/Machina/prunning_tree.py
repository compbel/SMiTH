"""
Author: Maryam KafiKang
Created: 2024
Description: This script processes phylogenetic trees for MACHINA analysis by pruning and
            reformatting them. It performs several key operations:
            1. Tree Pruning: Collapses nodes with identical host IDs to simplify the tree
            2. Node Renaming: Extracts and uses only host IDs from leaf labels
            3. Format Conversion: Generates MACHINA-compatible input files:
               - patientTree.txt: Edge list with unique IDs for nodes
               - patientLabeling.txt: Mapping of leaf IDs to original names
               - coloring.txt: Numerical mapping for unique leaf names

Input: Newick format phylogenetic trees from FAVITES output
Output: Pruned trees and MACHINA-compatible input files

Usage:
    The script processes multiple trees in batch, reading from:
    'log_50/FAVITES_output_logSI_contemp_T200_N100_E1_X/error_free_files/phylogenetic_trees/Tnet_Tree.0'
    where X ranges from 1 to 500.
"""

import os
import matplotlib.pyplot as plt
from Bio import Phylo

class PhylogenyPruner:
    def __init__(self, newick_file, output_filename, output_folder="pruned_trees_log50"):
        """
        Initializes the PhylogenyPruner class.

        Parameters:
        - newick_file (str): Path to the input Newick file.
        - output_folder (str): Directory where the pruned tree will be saved.
        - output_filename (str): Name of the pruned tree file.
        """
        self.newick_file = newick_file
        self.output_folder = f"{output_folder}/{output_filename}"
        self.output_filename = output_filename
        self.tree = self.load_tree()

    def load_tree(self):
        """Loads the phylogenetic tree from the Newick file."""
        return Phylo.read(self.newick_file, "newick")

    def extract_host_id(self, label):
        """
        Extracts the HostID from a label formatted as 'HostID_seqID'.
        If the label is missing or improperly formatted, return None.
        """
        if label and "_" in label:
            return label.split("_")[0]  # Extract HostID
        return None  # Return None for nodes without a valid HostID

    def rename_leaves(self, clade):
        """Renames all leaf nodes to only contain the HostID (removes seqID)."""
        if not clade.clades:  # If it's a leaf node
            clade.name = self.extract_host_id(clade.name)  # Update to HostID only
        else:
            for child in clade.clades:
                self.rename_leaves(child)

    def prune_tree(self, clade):
        """
        Recursively prunes the tree based on HostID.
        If all children have the same HostID (and not None), collapse them into the parent node.
        """
        if not clade.clades:  # If it's a leaf, return its HostID
            return clade.name

        child_host_ids = set(self.prune_tree(child) for child in clade.clades)

        # Remove None values from the set
        filtered_host_ids = {host_id for host_id in child_host_ids if host_id is not None}

        # If all children have the same non-None HostID, merge them
        if len(filtered_host_ids) == 1 and None not in child_host_ids:
            clade.name = list(filtered_host_ids)[0]  # Assign HostID to parent
            clade.clades = []  # Remove children
        return clade.name

    def plot_tree(self, title):
        """Plots the phylogenetic tree using Matplotlib."""
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(1, 1, 1)
        Phylo.draw(self.tree, axes=ax, do_show=False)
        plt.title(title)
        plt.show()

    def save_pruned_tree(self):
        """
        Saves the pruned tree in Newick format into the specified folder.
        Creates the folder if it does not exist.
        """
        os.makedirs(self.output_folder, exist_ok=True)
        output_path = os.path.join(self.output_folder, f"{self.output_filename}.nwk")
        Phylo.write(self.tree, output_path, "newick")
        #print(f"Pruned tree saved to: {output_path}")

    def generate_machina_input(self):
        """
        Generates input files for the Machina model:
        - patientTree.txt: Stores edges (assigns unique IDs to internal nodes and leaves).
        - patientLabeling.txt: Stores leaf IDs with original names.
        - coloring.txt: Maps each unique leaf name to a unique number.
        """
        os.makedirs(self.output_folder, exist_ok=True)

        patient_tree_path = os.path.join(self.output_folder, "patientTree.txt")
        patient_labeling_path = os.path.join(self.output_folder, "patientLabeling.txt")
        coloring_path = os.path.join(self.output_folder, "coloring.txt")

        leaf_nodes = []    # List of tuples (leaf_id, original_name)
        edges = []
        leaf_count = 0
        internal_count = 0
        node_index = {}  # Map nodes to unique identifiers

        def traverse_tree(clade, parent_id=None):
            nonlocal leaf_count, internal_count

            if not clade.clades:  # Leaf node
                node_id = f"leaf_{leaf_count}"
                leaf_nodes.append((node_id, clade.name))
                node_index[clade.name] = node_id
                leaf_count += 1
            else:  # Internal node
                node_id = f"internal_{internal_count}"
                internal_count += 1
                node_index[clade.name] = node_id

            if parent_id is not None:
                edges.append((parent_id, node_id))

            for child in clade.clades:
                traverse_tree(child, node_id)

        traverse_tree(self.tree.clade)

        # Write patientTree.txt (edges)
        with open(patient_tree_path, "w") as f:
            for parent, child in edges:
                f.write(f"{parent} {child}\n")

        # Write patientLabeling.txt (leaf ID with original name)
        with open(patient_labeling_path, "w") as f:
            for node_id, original_name in leaf_nodes:
                f.write(f"{node_id} {original_name}\n")

        # Create a dictionary to assign a unique number to each unique leaf name
        unique_leaf_numbers = {}
        for _, original_name in leaf_nodes:
            if original_name not in unique_leaf_numbers:
                unique_leaf_numbers[original_name] = len(unique_leaf_numbers)

        # Write coloring.txt with unique leaf names and their assigned number
        with open(coloring_path, "w") as f:
            for leaf_name, number in unique_leaf_numbers.items():
                f.write(f"{leaf_name} {number}\n")

        #print(f"Machina input files saved in {self.output_folder}")

    def process_tree(self):
        """Executes the full pipeline: rename, prune, save, and generate Machina input."""
        self.rename_leaves(self.tree.clade)
        self.prune_tree(self.tree.clade)
        self.save_pruned_tree()
        self.generate_machina_input()


for i in range(1,501):
    newick_file_path = f"log_50/FAVITES_output_logSI_contemp_T200_N100_E1_{i}/error_free_files/phylogenetic_trees/Tnet_Tree.0"
    # Create an instance of PhylogenyPruner
    if os.path.exists(newick_file_path):
        pruner = PhylogenyPruner(newick_file_path, output_filename=f"FAVITES_output_logSI_contemp_T200_N100_E1_{i}")
        # Run the entire process
        pruner.process_tree()