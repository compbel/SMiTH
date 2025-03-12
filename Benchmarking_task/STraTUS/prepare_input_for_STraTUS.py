import os
import csv
from Bio import Phylo

def extract_host_tip_mapping(newick_file, output_csv):
    """
    Reads a phylogenetic tree in Newick format, extracts tip labels,
    and generates a CSV file mapping each tip to its host.

    Parameters:
    - newick_file: Path to the input phylogenetic tree in Newick format.
    - output_csv: Path to the output CSV file for STraTUS.
    """

    # Load the phylogenetic tree
    tree = Phylo.read(newick_file, "newick")

    # Extract tip labels and map to host ID
    tip_host_mapping = []
    for tip in tree.get_terminals():
        if "_" in tip.name:
            host_id, seq_id = tip.name.split("_", 1)  # Extract Host ID from format: HostID_SeqID
            tip_host_mapping.append((tip.name, f'{host_id}H'))  # Store (Tip Label, Host ID)

    # Write to CSV file
    with open(output_csv, mode="w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["Tip Label", "Host ID"])  # Header
        writer.writerows(tip_host_mapping)  # Write data

    #print(f"CSV file created: {output_csv}")

# Example usage
for i in range(1,501):
    newick_tree_file = f"Favites_log/FAVITES_output_logSI_contemp_T200_N100_E1_{i}/error_free_files/phylogenetic_trees/Tnet_Tree.0"  # Replace with your actual Newick file
    output_csv_file = f"host_data/FAVITES_output_logSI_contemp_T200_N100_E1_{i}.csv"
    if os.path.exists(newick_tree_file):
        extract_host_tip_mapping(newick_tree_file, output_csv_file)
