"""
Author: Maryam KafiKang
Description: This script evaluates the performance of the Phyloscanner tool by comparing its inferred
            transmission networks against ground truth data from FAVITES. It calculates key metrics
            such as sensitivity, specificity, and F1 score, and records execution times for analysis.
            The workflow includes:
            1. Data Reading:
               - Loads ground truth and inferred networks
            2. Metric Calculation:
               - Computes sensitivity, specificity, and F1 score
            3. Output Generation:
               - Saves detailed statistics and execution times to files

Input:
    - Ground truth networks from FAVITES
    - Inferred networks from Phyloscanner

Output:
    - CSV files with performance metrics
    - Execution time statistics
"""

import os
import sys
from os import listdir
from os.path import join, isfile

import re

metrics = ['sensitivity', 'specificity', 'f1']
tools = ['phyloscanner']



def read_ground_truth_network(folder):
    path = join(folder, 'error_free_files/transmission_network.txt')
    with open(path) as file:
        lines = file.readlines()
    edges = set()
    for line in lines[1:]:  # Skip the root node
        u, v = map(int, line.split('\t')[:2])
        if u != v:
            edges.add((u, v))
    return edges


def read_network(folder):
    edges = set()
    running_time = 0.0
    with open(folder, 'r') as file:
        for line in file:
            line = line.strip()
            if not line or "Running Time" in line:
                if "Running Time" in line:
                    running_time = float(re.search(r'\d+\.\d+', line).group())
                continue
            sp = line.split('\t')
            if len(sp) >= 2:
                edge = (int(sp[0]), int(sp[1]))
                edges.add(edge)
    return edges, running_time


def calculate_stat(ground_truth, inferred):
    true_positives = len(ground_truth & inferred)
    sensitivity = true_positives / len(ground_truth) if ground_truth else 0
    specificity = true_positives / len(inferred) if inferred else 0
    f1 = (2 * sensitivity * specificity / (sensitivity + specificity)) if (sensitivity + specificity) > 0 else 0
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
    out = open('phyloscanner_results.txt', 'w')
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

    output_file = "phyloscanner_execution_time.txt"

    # Write dictionary to file
    with open(output_file, 'w') as file:
        for key, value in time.items():
            file.write(f"{key}: {value:.5f}\n")
    

categories = {}
time = {}
# reads the outputs for all tools and calculate metrics with regard to the ground truth
def run_analysis(folder, tnet_folder):
    
    sample_name = folder.split('/')[-2]
    config = folder.split('/')[-3]
    
    category = config
    ground_truth = read_ground_truth_network(folder)
    tnet, running_time = read_network(tnet_folder)
    time[sample_name] =  running_time
    

    if len(ground_truth) != 0:
        if category not in categories:
            categories[category] = {metric: [] for metric in metrics}
            categories[category]['sample_names'] = []
            categories[category]['time'] = []

        values = calculate_stat(ground_truth, tnet)
        for metric, value in zip(metrics, values):
            categories[category][metric].append(value)
        categories[category]['sample_names'].append(sample_name)
        categories[category]['time'].append(running_time)

# example to run python3 favites_validation.py Favites_inputs/ data
if __name__ == "__main__":
    base_path = "Favites_log"
    for i in range(1, 501):
        print(i)
        true_folder = f"{base_path}/FAVITES_output_logSI_contemp_T200_N100_E1_{i}/"
        generated_folder = f"phyloscanner_output/FAVITES_output_logSI_contemp_T200_N100_E1_{i}_phyloscannerTree.txt"
        if os.path.exists(generated_folder):
            run_analysis(true_folder, generated_folder)
    print_results(categories, time)

