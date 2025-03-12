import csv
import os
import pandas as pd
import subprocess
import generate_TransmissionTree
import time
import argparse
from datetime import datetime

class Logger(object):
    def __init__(self, path, header, write_header=True):
        """
        Create an instance of the logger.
        Parameters
        ----------
        path: str
            The file path where the log file will be written.
        header: list
            The header line of the log file.
        """
        self.log_file = open(path, 'w')  # Use 'wb' for Python 2 compatibility
        self.logger = csv.writer(self.log_file, delimiter='\t')

        if write_header:
            self.logger.writerow(header)
        self.header = header

    def __del__(self):  # Fixed method name
        self.log_file.close()

def run_phyloscanner(i, file_name, output_dir, modes=[0]):
    command = """Rscript phyloscanner_analyse_trees.R {} {} s,0 -x "(.*)_N([0-9]+)$" --allClassifications -od {}""".format(file_name, i, output_dir)
    #command = """Rscript phyloscanner_analyse_trees.R {0} {1} s,0 -x "(.*)_N([0-9]+)$" --allClassifications """.format(file_name, output_dir)
    #print(command)
    os.system(command)

if __name__ == '__main__':
    #### save execution time
    parser = argparse.ArgumentParser(description="Run Phyloscanner on FAVITES output")
    parser.add_argument("i", type=int, help="Index parameter")
    args = parser.parse_args()

    timing_dir = "phyloscanner_execution_time.txt"
    train_logger = Logger(
        timing_dir,
        ['Test Data', 'Model', 'Execution Time (seconds)', '#hosts']
    )

    data = "log_50"
    
    #for i in range(1, 2):
    i = args.i
    path = data + "/FAVITES_output_logSI_contemp_T200_N100_E1_" + str(i) + "/error_free_files/phylogenetic_trees/Tnet_Tree.0"
    
    if os.path.isfile(path):
        dir = "phyloscanner_output_phylo/"
        if not os.path.exists(dir):
            os.makedirs(dir)  # Removed exist_ok=True (Python 3 feature)

        output_dir = "phyloscanner_output/"
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        start_time = time.time()
        run_phyloscanner("FAVITES_output_logSI_contemp_T200_N100_E1_" + str(i), path, output_dir=dir, modes=[0])
        
        
        if os.path.exists(dir + "FAVITES_output_logSI_contemp_T200_N100_E1_" + str(i) + "_classification.csv"):
            ts_tree_file = os.path.join(output_dir, "FAVITES_output_logSI_contemp_T200_N100_E1_" + str(i) + "_phyloscannerTree.txt")
            hosts = generate_TransmissionTree.create_transmission_network(dir + "FAVITES_output_logSI_contemp_T200_N100_E1_" + str(i) + "_classification.csv", ts_tree_file)
            end_time = time.time()
            execution_time = end_time - start_time
            with open(ts_tree_file, 'a') as f:
                f.write("\nRunning Time (seconds): {:.5f}\n".format(execution_time))
