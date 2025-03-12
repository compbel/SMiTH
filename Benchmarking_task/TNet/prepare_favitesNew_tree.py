"""
Author: Maryam KafiKang
Created: 2024
Description: This script prepares FAVITES phylogenetic tree files for TNet analysis
            by converting tree labels from FAVITES format to TNet compatible format.
            It processes multiple simulation outputs in batch.
"""

import os
import re
from pathlib import Path



def rename_favites_ph_tree_labels_into_tnet_labels(file_name):
    """
    Converts FAVITES tree labels from the format 'N58|3|100.0:' to '3_N58:'
    """
    output_file_name = file_name+"Tnet_Tree.0"

    with open(file_name +'tree_0.time.tre', 'r') as f, open(output_file_name, 'w') as outf:
        data = re.sub(r'(N\d+)\|(\d+)\|(\d+\.\d+):', r'\2_\1:', f.read())
        outf.write(data)
    

def read_files_from_directories(root_dir):
    
    # Make sure that path is a Path object
    path = f'{root_dir}/FAVITES_output_logSI_contemp_T200_N100_E1_'  # Using pathlib to handle path operations correctly
    #if path.is_dir():  # Checking if path is indeed a directory

    for i in range(1,501):
        #input_file =  path + '/phylogeny_contracted_newick_' + str(i)  
        
        file_path =  path + str(i) 
        # Assuming you want to append the specific path to file_path
        input_file = f'{file_path}/error_free_files/phylogenetic_trees/'
        if os.path.exists(f"{input_file}/tree_0.time.tre"):
            if os.path.isdir(file_path):
                rename_favites_ph_tree_labels_into_tnet_labels(input_file)
        else:
            print(f'no tree for {file_path}')

root_directory = 'log_50'

read_files_from_directories(root_directory)

