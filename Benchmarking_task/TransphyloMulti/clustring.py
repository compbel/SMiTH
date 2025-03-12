"""
Author: Maryam KafiKang
Description: This script performs sequence clustering and tree pruning for phylogenetic analysis.
            It processes sequence data to identify representative sequences for each host and prunes
            phylogenetic trees accordingly. The workflow includes:
            1. Sequence Parsing:
               - Reads sequences from FASTA files
               - Organizes sequences by host
            2. Feature Extraction:
               - Computes k-mer frequency features for each sequence
            3. Clustering:
               - Clusters sequences for each host using k-means
               - Identifies representative sequences for each cluster
            4. Tree Pruning:
               - Prunes phylogenetic trees to retain only representative sequences
Input:
    - FASTA files with sequence data
    - Newick format phylogenetic trees

Output:
    - Clustered representative sequences in FASTA format
    - Pruned phylogenetic trees in Newick format
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import numpy as np
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.cluster import KMeans
from scipy.spatial.distance import cdist
from ete3 import Tree
import os

class SequenceClustering:
    def __init__(self, fasta_file, output, n_clusters=10, k_mer=3):
        """
        Initialize the SequenceClustering class.
        
        Args:
            fasta_file (str): Path to the input FASTA file
            n_clusters (int): Number of clusters per host
            k_mer (int): Length of k-mers for feature extraction
        """
        self.fasta_file = fasta_file
        self.n_clusters = n_clusters
        self.k_mer = k_mer
        self.n_init = 10
        self.output_file = f"sequence{output}.txt"
        self.sequence_info = defaultdict(dict)  # Store sequence info by host
        self.cluster_seq = []  # Store cluster info by host

    def parse_fasta(self):
        """
        Parse the FASTA file and organize sequences by host.
        Returns a nested structure:
        {
            host_id: {
                'sequences': [...],  # actual sequences
                'seq_ids': [...],    # sequence IDs
                'features': None     # will store k-mer features later
            }
        }
        """
        for record in SeqIO.parse(self.fasta_file, "fasta"):
            # Parse the header: >N1|1|100
            seq_id, host_id, _ = record.id.split('|')
            
            if host_id not in self.sequence_info:
                self.sequence_info[host_id] = {
                    'sequences': [],
                    'seq_ids': [],
                    'features': None
                }
            
            self.sequence_info[host_id]['sequences'].append(str(record.seq))
            self.sequence_info[host_id]['seq_ids'].append(seq_id)
        
        # Convert lists to numpy arrays
        for host_data in self.sequence_info.values():
            host_data['sequences'] = np.array(host_data['sequences'])
            host_data['seq_ids'] = np.array(host_data['seq_ids'])

    def compute_features(self):
        """
        Compute k-mer frequency features for all sequences of each host.
        """
        vectorizer = CountVectorizer(analyzer='char', ngram_range=(self.k_mer, self.k_mer))
        
        for host_data in self.sequence_info.values():
            features = vectorizer.fit_transform(host_data['sequences'])
            host_data['features'] = features.toarray()

    def cluster_sequences(self):
        """
        Cluster sequences for each host and store representative sequences.
        Updates sequence_info with cluster assignments and representatives.
        """
        for host_id, host_data in self.sequence_info.items():
            # Perform clustering
            kmeans = KMeans(n_clusters=self.n_clusters, random_state=42, n_init=self.n_init)
            labels = kmeans.fit_predict(host_data['features'])
            
            # Initialize cluster information
            host_data['clusters'] = defaultdict(list)
            host_data['representatives'] = defaultdict(dict)
            
            # Compute centroids
            centroids = np.vstack([host_data['features'][labels == i].mean(axis=0) 
                                 for i in range(self.n_clusters)])
            
            # For each cluster, find representative sequence
            for cluster_id in range(self.n_clusters):
                cluster_mask = (labels == cluster_id)
                cluster_features = host_data['features'][cluster_mask]
                cluster_sequences = host_data['sequences'][cluster_mask]
                cluster_seq_ids = host_data['seq_ids'][cluster_mask]
                
                # Store all sequence IDs in this cluster
                host_data['clusters'][cluster_id] = cluster_seq_ids.tolist()
                
                # Find sequence closest to centroid
                distances = cdist([centroids[cluster_id]], cluster_features, metric="euclidean").flatten()
                rep_idx = np.argmin(distances)
                
                # Store representative information
                host_data['representatives'][cluster_id] = {
                    'seq_id': cluster_seq_ids[rep_idx],
                    'sequence': cluster_sequences[rep_idx],
                    'cluster_size': np.sum(cluster_mask)
                }

    def write_output(self, output_file):
        """
        Write representative sequences to a FASTA file and cluster information to a summary file.
        """
        # Write FASTA file with representatives
        records = []
        for host_id, host_data in self.sequence_info.items():
            for cluster_id, rep_info in host_data['representatives'].items():
                seq_id = f"{host_id}_{rep_info['seq_id']}"
                self.cluster_seq.append(seq_id)
                record = SeqRecord(Seq(rep_info['sequence']), id=seq_id, description="")
                records.append(record)

        # Create output directory if it doesn't exist
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        
        # Write representative sequences
        with open(output_file, "w") as output_handle:
            SeqIO.write(records, output_handle, "fasta")
        return self.cluster_seq
        
        # # Write cluster summary
        # summary_file = self.output_file.replace('.fasta', '_summary.txt')
        # with open(summary_file, 'w') as f:
        #     for host_id, host_data in self.sequence_info.items():
        #         f.write(f"\nHost {host_id} clusters:\n")
        #         for cluster_id, seq_ids in host_data['clusters'].items():
        #             rep_info = host_data['representatives'][cluster_id]
        #             f.write(f"\nCluster {cluster_id}:\n")
        #             f.write(f"Representative: {rep_info['seq_id']}\n")
        #             f.write(f"Cluster size: {rep_info['cluster_size']}\n")
        #             f.write(f"All sequences: {', '.join(seq_ids)}\n")

    def main(self):
        """
        Main workflow for sequence clustering and representative selection.
        """
        self.parse_fasta()
        self.compute_features()
        self.cluster_sequences()
        
        # Generate output filename
        output_dir = "ClusteredData"
        #basename = os.path.basename(self.fasta_file)
        output_file = os.path.join(output_dir, self.output_file)
        
        cluster_dic = self.write_output(output_file)
        return cluster_dic

# Example usage
if __name__ == "__main__":
    directory = "Favites_exp"

    # Dictionary to organize sequences by host
    sequences_by_host = {}
    folders = [
        folder
        for folder in os.listdir(directory)
        if os.path.isdir(os.path.join(directory, folder)) and folder.startswith("FAVITES_output")
    ]
    for folder in folders:
        fasta_file = f'{directory}/{folder}/error_free_files/sequence_data.fasta'
        output_num = folder.split("N100")[1]
        clustering = SequenceClustering(fasta_file, output_num, n_clusters=10, k_mer=3)
        representative_seq = clustering.main()
        #time for prunning rest of tree and keep cluster_dic
        tree = Tree(f"{directory}/{folder}/error_free_files/phylogenetic_trees/Tnet_Tree.0")

        # List of leaf names you want to keep
        leaves_to_keep = representative_seq  # Replace with your sequence IDs

        # Prune the tree
        tree.prune(leaves_to_keep, preserve_branch_length=True)

        outfile=f"tree/PrunedTree{output_num}.nwk"
        # Write the pruned tree to file
        os.makedirs(os.path.dirname(outfile), exist_ok=True)
        tree.write(format=1, outfile = outfile)
    
