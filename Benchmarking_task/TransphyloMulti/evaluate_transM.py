"""
Author: Maryam KafiKang
Description: This script evaluates transmission network inference models by comparing generated
            networks against ground truth data. It calculates various performance metrics based on the different thresholds using precision-recall curves
            and generates visualizations to assess model accuracy. The workflow includes:
            1. Data Reading:
               - Loads ground truth and inferred networks
               - Reads intermediate network scores for detailed analysis
            2. Metric Calculation:
               - Computes sensitivity, specificity, and F1 score
               - Evaluates performance at multiple probability thresholds
            3. Visualization:
               - Generates precision-recall curves for individual samples
               - Saves plots for comparative analysis

Input:
    - Ground truth networks from FAVITES
    - Inferred networks from various models
    - Intermediate network scores for detailed evaluation

Output:
    - CSV files with performance metrics
    - Precision-recall curve plots
    - Execution time statistics
"""

import os
from os.path import join, isfile
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc, precision_recall_curve
import seaborn as sns
import math
import numpy as np
from collections import defaultdict
import pandas as pd

class NetworkAnalysis:
    def __init__(self, tool, True_network_dir, generated_tranmission_network_dir, genearetd_intermediate_network_dir):
        self.tool = tool
        self.True_network_dir = True_network_dir
        self.generated_tranmission_network_dir = generated_tranmission_network_dir
        self.genearetd_intermediate_network_dir = genearetd_intermediate_network_dir
        self.metrics = ['sensitivity', 'specificity', 'f1']
        self.categories = {}
        self.time = {}
        self.sensitivity, self.specificity, self.f1= [],[],[]
        self.true_edge = []
        self.gen_edge = []
        self.sample_name = []
        self.threshlds = [0, 0.1, 0.2, 0.3, 0.4, 0.5]

    
    def read_ground_truth_network(self, folder):
        path = f'{folder}/error_free_files/transmission_network.txt'
        if not isfile(path):
            return []
        lines = open(path).readlines()
        edges = [(int(sp[0]), int(sp[1])) for line in lines[1:] if (sp := line.split('\t')) and int(sp[0]) != int(sp[1])]
        return edges
    
    def read_intermediate_network(self, sample_name, filepath):
        if not isfile(filepath):
            return []

        with open(filepath, "r") as file:
            lines = file.readlines()

        # Extract the column headers from the first line
        headers = list(map(int, lines[0].strip().split()))

        # Prepare the data rows (skip first row since it's just column headers)
        data = []
        for line in lines[1:]:
            values = list(map(float, line.strip().split()))
            row_label = int(values[0])  # First value is the row index
            row_data = values[1:]  # Remaining values are the row contents
            data.append([row_label] + row_data)

        # Convert to DataFrame
        df = pd.DataFrame(data, columns=["ID"] + headers)
        df.set_index("ID", inplace=True)  # Set row names


        return df


    def read_network(self,sample_name, filepath):
        if not isfile(filepath):
            return []

        with open(filepath, "r") as file:
            lines = file.readlines()

        if "Running time" in lines[-1]:
            time = "".join([c for c in lines[-1] if c.isdigit() or c == "."])

        # Convert list of lines back to a string for pandas
        cleaned_data = "".join(lines[:-1])

        # Read the cleaned data into a DataFrame
        df = pd.read_csv(pd.io.common.StringIO(cleaned_data), sep="\t", index_col=0)

        #now we have to convert the dataframe to a list of tuples in a way that has this info: row, column, max probability, average probability
        transmission_links = []
        for infector in df.index:
            for infectee in df.columns:
                prob = df.loc[infector, infectee]
                if prob > 0:  # Only consider probabilities greater than zero
                    transmission_links.append((infector, int(infectee), prob))

        result = []
        for index, row in df.iterrows():
            non_zero_values = row[row!=0]
            if len(non_zero_values) != 0:
                max_prob = non_zero_values.max()
                avg_prob = non_zero_values.mean()
            else:
                max_prob = 0
                avg_prob = 0
            col = row.idxmax()
            result.append((index, col, max_prob, avg_prob))


        base_name = os.path.basename(self.generated_tranmission_network_dir)

        # Construct the output directory
        output_dir = os.path.join(
            os.path.dirname(self.generated_tranmission_network_dir),  # Parent directory
            f"{base_name}_Extracted"  # Append "_Extracted" to the directory name
        )

        # Construct the full output file path
        output_path = os.path.join(output_dir, f"{sample_name}_extracted_network.txt")

        # Ensure the directory exists
        os.makedirs(output_dir, exist_ok=True)

        #write the result info as transmission network 
        with open(output_path, 'w') as out:
            
            out.write("First Header: Transmission Data\n")
            out.write("Infector\tInfectee\tProbability\n")  # Column names
            for row, max_val, avg_val in transmission_links:
                out.write(f"{row}\t{max_val}\t{avg_val}\n")
            out.write("----------------------------------------------\n")  # Separator
            
            # Second header
            out.write("Second Header: Summary Statistics\n")
            out.write("Infector\tHighProbableInfectee\tMaxProbabily\tAvg Probability (P>0)\n")  # Second section column names
            for row, col, max_val, avg_val in result:
                out.write(f"{row}\t{col}\t{max_val}\t{avg_val}\n")
            
        
        return transmission_links, time

    def calculate_stat(self, ground_truth, inferred):
        sensitivity = len([e for e in ground_truth if e in inferred]) / len(ground_truth) if ground_truth else 0
        specificity = len([e for e in inferred if e in ground_truth]) / len(inferred) if inferred else 0
        f1 = (2 * sensitivity * specificity / (sensitivity + specificity)) if (sensitivity + specificity) > 0 else 0
        return sensitivity, specificity, f1

    def print_results(self):
        with open(f'{self.tool}_results.csv', 'w') as out:
            out.write("Sample name\t")
        
            headers = []
            for th in self.threshlds: 
                headers.append(f"Sensitivity_{th}")
                headers.append(f"Specificity_{th}")
                headers.append(f"F1_{th}")
                headers.append(f"intermediate_points_{th}")

            # Write the modified headers
            out.write("\t".join(headers) + "\tRunning_time(sec)" + "\n")

            num_expected_values = len(headers)  # The expected number of values for each sample

            for sample_name, data in self.categories.items():
                values = [sample_name]  # Start with the sample name
                
                for category, metrics in data.items():
                    for key, val in metrics.items():
                        values.append(val[0])  # Append the metric values
                
                # Fill missing values with zeros if not enough values are present
                while len(values) - 1 < num_expected_values:  # Subtracting 1 for sample_name
                    values.append(0)

                # Write the modified values
                out.write("\t".join(map(str, values)) + f'\t{self.time.get(sample_name, 0)}' + "\n")

    def plot_individual_pr_curves(self, true_edges_list, generated_edges_list, sample, max_cols=5):
        num_datasets = len(true_edges_list)

        
        cols = min(max_cols, num_datasets)  # Limit the number of columns (max 5 for readability)
        rows = math.ceil(num_datasets / cols)  # Ensure all datasets fit in grid

        fig, axes = plt.subplots(rows, cols, figsize=(cols * 4, rows * 4))  # Adjust figure size dynamically
        axes = np.array(axes).flatten()

        for i, (true_edges, generated_edges) in enumerate(zip(true_edges_list, generated_edges_list)):
            true_set = set(true_edges)
            gen_edge_dict = {(e1, e2): score for e1, e2, score in generated_edges}

            y_true = []
            y_pred = []
            
            for edge in true_set:
                y_true.append(1)
                y_pred.append(gen_edge_dict.get(edge, 0))
            
            for edge in gen_edge_dict:
                if edge not in true_set:
                    y_true.append(0)
                    y_pred.append(gen_edge_dict[edge])
            
            precision, recall, _ = precision_recall_curve(y_true, y_pred)
            pr_auc = auc(recall, precision)

            ax = axes[i]
            ax.plot(recall, precision, label=f'AUC = {pr_auc:.2f}', color='blue')
            ax.set_xlabel('Recall')
            ax.set_ylabel('Precision')
            ax.set_title(f'Sample {sample[i].split("_E1_")[1]}')
            ax.legend()
            ax.grid()


        for j in range(i + 1, len(axes)):
            fig.delaxes(axes[j])

        plt.tight_layout()
        plt.savefig('PR_curve.png', dpi=300)
        plt.show()

    def run_analysis(self, folder, tool_folder, intermediate_network):
        sample_name = folder.split('/')[-1]
        #category = folder.split('/')[-2]
        
        ground_truth = self.read_ground_truth_network(folder)
        inferred, time = self.read_network(sample_name, tool_folder)
        intermediate_df = self.read_intermediate_network(sample_name, intermediate_network)

        self.time[sample_name] = time
        
        
        infered_threshold = defaultdict(list)
        intermediate_score_dict = {}
        for t in self.threshlds:
            intermediate_score = 0
            for e1, e2, score in inferred:
                if score >= t:
                    infered_threshold[t].append((e1, e2))
                    intermediate_score += intermediate_df.loc[e1, e2] 
            intermediate_score_dict[t] = intermediate_score / len(inferred)


        if ground_truth:
            for category, inferred_pairs in infered_threshold.items():
                #for sample_name in inferred_pairs:  # Assuming inferred_pairs is sample-specific
                if sample_name not in self.categories:
                    self.categories[sample_name] = {}  # Make sample names the top-level keys

                if category not in self.categories[sample_name]:
                    self.categories[sample_name][category] = {metric: [] for metric in self.metrics}

                values = self.calculate_stat(ground_truth, inferred_pairs)  

                for metric, value in zip(self.metrics, values):
                    self.categories[sample_name][category][metric].append(round(value, 3)) 
                self.categories[sample_name][category]['intermediate_points'] = [intermediate_score_dict[category]]
            

    def write_running_time(self):
        with open(f'{self.tool}_running_time.txt', 'w') as out:
            for sample_name, time in self.time.items():
                out.write(f"{sample_name}: {time}\n")
        
    def process_configs(self):
        folders = [folder.split("_transphylo.txt")[0] for folder in os.listdir(self.generated_tranmission_network_dir) if folder != '.DS_Store']
        for folder in folders:
            generated_network_folder = f"{self.generated_tranmission_network_dir}/{folder}_transphylo.txt"
            true_network_folder = f"{self.True_network_dir}/{folder}"
            genearetd_intermediate_network = f"{self.genearetd_intermediate_network_dir}/{folder}_interMatrix.txt"
            
            print(f"{folder} is processing ...")
            self.run_analysis(true_network_folder, generated_network_folder, genearetd_intermediate_network)
        self.print_results()
        #self.plot_individual_pr_curves(self.true_edge, self.gen_edge, self.sample_name)
        self.write_running_time()

if __name__ == "__main__":

    generated_tranmission_network_dir = "TransphyloMultiResultMat/Shape_0.1_Scale_200/wiwMatrix" #Shape_0.1_Scale_200
    genearetd_intermediate_network_dir = "TransphyloMultiResultMat/Shape_0.1_Scale_200/intermediateMatrix"
    Truth_network_dir = "Favites_exp" #Data
    tool_name = "TransphyloMulti_FavitesExp"  # Change this to use another tool
    analysis = NetworkAnalysis(tool_name, Truth_network_dir, generated_tranmission_network_dir,genearetd_intermediate_network_dir )
    analysis.process_configs()
    

    
