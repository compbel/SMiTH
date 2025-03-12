import os
import sys
from os import listdir
from os.path import join, isfile
import numpy as np

import re

metrics = ['sensitivity', 'specificity', 'f1']
tools = ['machina']



def read_ground_truth_network(folder):
    path = join(folder, 'error_free_files/transmission_network.txt')
    lines = open(path).readlines()
    edges = []
    sp = lines[0].split('\t')
    # skip the root since phyloscanner doesn't give the root
    # edges.append((-1, int(sp[1])))
    for line in lines[1:]:
        sp = line.split('\t')
        u = int(sp[0])
        v = int(sp[1])
        if u == v:
            continue
        edges.append((u, v))
    return edges


def read_network(folder):
    lines = open(folder).readlines()
    edges = []
    for line in lines[0:]:
        sp = line.split('\t')[0].split('\n')[0].split(" ")
        edge = (int(sp[0]), int(sp[1]))
        if edge not in edges:
            edges.append(edge)
    return edges


# calculates the metrics for reconstructed networks by comparing given edge lists
def calculate_stat(ground_truth, inferred):
    sensitivity = len([e for e in ground_truth if e in inferred]) / len(ground_truth)
    if len(inferred) == 0:
        specificity = 0
    else:
        specificity = len([e for e in inferred if e in ground_truth]) / len(inferred)
    if sensitivity < 0.00001 or specificity < 0.00001:
        f1 = 0
    else:
        f1 = 2 * (sensitivity * specificity) / (sensitivity + specificity)
    return sensitivity, specificity, f1

# prints statistics to console and save detailed statistics by categories in results.txt
def print_results(results, time):
    '''''
    for category in results:
        n_samples = len(results[category]['sensitivity'])
        print("Category {}. Total samples {}".format(category, n_samples))

        for tool in tools:
            values = [sum(results[category][metric]) / n_samples for metric in metrics]
            print(tool + " ", end='')
            print(",".join(["{} = {:.3f}".format(metric, value) for (metric, value) in zip(metrics, values)]))
    '''''
    out = open(f"{tools[0]}_results.txt", 'w')
    for category in results:
        out.write("category:{}\n".format(category))
        column_names = ['sample_name', 'sensitivity', 'specificity', 'f1']

        out.write(",".join(column_names))
        out.write("\n")
        sample_n = 0
        for sample_name in results[category]['sample_names']:
            values = [sample_name]
            for tool in tools:
                for metric in metrics:
                    values.append("{:.3f}".format(results[category][metric][sample_n]))
            sample_n += 1
            out.write(",".join(values))
            out.write("\n")

    with open(f"{tools[0]}_runing_times.txt", 'w') as file:
        for key, value in time.items():
            file.write(f"{key}: {value}\n")



categories = {}
# reads the outputs for all tools and calculate metrics with regard to the ground truth
def run_analysis(folder, model_folder):
    
    sample_name = folder.split('/')[-2]
    config = folder.split('/')[-3]
    
    category = config
    ground_truth = read_ground_truth_network(folder)
    model = read_network(model_folder)
    

    if len(ground_truth) != 0:
        if category not in categories:
            categories[category] = {metric: [] for metric in metrics}
            categories[category]['sample_names'] = []

        values = calculate_stat(ground_truth, model)
        for metric, value in zip(metrics, values):
            categories[category][metric].append(value)
        categories[category]['sample_names'].append(sample_name)


def extract_execution_time(filename):
    """
    Reads a text file and extracts the execution time from the second line
    that contains 'PS, S'.
    """
    execution_time = None  # Variable to store execution time
    pattern = re.compile(r"PS, S")  # Regex pattern to match 'PS, S'
    
    with open(filename, "r") as file:
        for line in file:
            if pattern.search(line):  # Find the line that contains 'PS, S'
                parts = line.strip().split()  # Split the line by spaces/tabs
                execution_time = float(parts[-1])  # Get the last number
                break  # Stop after finding the first match

    return execution_time

execution_time = {}
if __name__ == "__main__":
    pattern = re.compile(r"^\d+-G-\d+-S\.tree$")  # Regex pattern to match 'number_G_number_S.txt'

    for i in range(1, 501):
        fav_folder = f"Favites_log/FAVITES_output_logSI_contemp_T200_N100_E1_{i}/"
        
        # List all files in the machina results folder
        model_folder_path = f"machina_results/FAVITES_output_logSI_contemp_T200_N100_E1_{i}/"
        if os.path.exists(model_folder_path):
            # Process only files that match the pattern
            for filename in os.listdir(model_folder_path):
                if pattern.match(filename):  # Check if filename matches 'number_G_number_S.txt'
                    model_folder = os.path.join(model_folder_path, filename)  # Full file path
                    
                    if os.path.isfile(model_folder):  # Ensure it's a file
                        run_analysis(fav_folder, model_folder)
            time = extract_execution_time(f"{model_folder_path}result.txt")
            execution_time[f'FAVITES_output_logSI_contemp_T200_N100_E1_{i}'] = time

    print_results(categories, execution_time)
