"""
Author: Maryam KafiKang
Description: This script evaluates STraTUS's performance by comparing inferred transmission networks
            with FAVITES ground truth networks. It calculates and reports various metrics including:
            - Sensitivity (true positive rate)
            - Specificity (true negative rate)
            - F1 score (harmonic mean of precision and recall)
"""

import os
import sys
from os.path import join
from typing import Set, Tuple, List, Dict, Optional
from collections import defaultdict
import re
import csv
from datetime import datetime
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# Type aliases for better code readability
Edge = Tuple[int, int]
EdgeSet = Set[Edge]
MetricResult = Tuple[float, float, float]

class StratusEvaluator:
    def __init__(self):
        self.metrics = ['sensitivity', 'specificity', 'f1']
        self.categories: Dict[str, Dict] = defaultdict(lambda: {
            'sensitivity': [], 'specificity': [], 'f1': [],
            'sample_names': [], 'time': []
        })

    @staticmethod
    def read_ground_truth_network(folder: str) -> EdgeSet:
        """Read and parse ground truth network from FAVITES output."""
        try:
            path = join(folder, 'error_free_files/transmission_network.txt')
            edges: EdgeSet = set()
            
            with open(path) as f:
                # Skip first line (root node)
                next(f)
                for line in f:
                    try:
                        u, v = map(int, line.strip().split('\t')[:2])
                        if u != v:  # Skip self-loops
                            edges.add((u, v))
                    except (ValueError, IndexError):
                        continue
            return edges
        except FileNotFoundError:
            logging.warning(f"Ground truth file not found in {folder}")
            return set()
        except Exception as e:
            logging.error(f"Error reading ground truth file: {e}")
            return set()

    @staticmethod
    def read_network(filepath: str) -> EdgeSet:
        """Read and parse inferred network file from STraTUS output."""
        edges: EdgeSet = set()
        
        try:
            with open(filepath) as file:
                for line in file:
                    line = line.strip()
                    if not line:
                        continue

                    parts = line.split()
                    if len(parts) < 3:  # Need at least source, target, and weight
                        continue

                    try:
                        # Only add edge if weight is non-zero
                        if float(parts[2]) != 0.0:
                            edge = (int(parts[0]), int(parts[1]))
                            edges.add(edge)
                    except (ValueError, IndexError):
                        continue

            return edges
        except FileNotFoundError:
            logging.warning(f"Network file not found: {filepath}")
            return set()
        except Exception as e:
            logging.error(f"Error reading network file {filepath}: {e}")
            return set()

    @staticmethod
    def calculate_metrics(ground_truth: EdgeSet, inferred: EdgeSet) -> MetricResult:
        """Calculate performance metrics using set operations for efficiency."""
        if not ground_truth:
            return 0.0, 0.0, 0.0

        true_positives = len(ground_truth.intersection(inferred))
        sensitivity = true_positives / len(ground_truth) if ground_truth else 0.0
        specificity = true_positives / len(inferred) if inferred else 0.0

        # Calculate F1 score
        if sensitivity < 0.00001 or specificity < 0.00001:
            f1 = 0.0
        else:
            f1 = 2 * (sensitivity * specificity) / (sensitivity + specificity)

        return sensitivity, specificity, f1

    def save_results(self, output_dir: str = '.') -> None:
        """Save results to CSV files for better data handling."""
        os.makedirs(output_dir, exist_ok=True)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

        # Save metrics
        metrics_file = join(output_dir, f'stratus_results_{timestamp}.csv')
        with open(metrics_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['category', 'sample_name', 'sensitivity', 'specificity', 'f1'])
            
            for category, data in self.categories.items():
                for idx, sample_name in enumerate(data['sample_names']):
                    writer.writerow([
                        category,
                        sample_name,
                        f"{data['sensitivity'][idx]:.3f}",
                        f"{data['specificity'][idx]:.3f}",
                        f"{data['f1'][idx]:.3f}"
                    ])

        logging.info(f"Results saved to {metrics_file}")

    def process_sample(self, true_folder: str, inferred_folder: str) -> None:
        """Process a single sample and calculate its metrics."""
        try:
            sample_name = os.path.basename(true_folder.rstrip('/'))
            category = os.path.basename(os.path.dirname(true_folder))

            ground_truth = self.read_ground_truth_network(true_folder)
            if not ground_truth:
                logging.warning(f"No ground truth edges found for {sample_name}")
                return

            inferred_edges = self.read_network(inferred_folder)
            if not inferred_edges:
                logging.warning(f"No inferred edges found for {sample_name}")
                return

            metrics = self.calculate_metrics(ground_truth, inferred_edges)
            cat_data = self.categories[category]
            
            for metric_name, value in zip(self.metrics, metrics):
                cat_data[metric_name].append(value)
            cat_data['sample_names'].append(sample_name)

            logging.info(f"Processed sample {sample_name} successfully")
        except Exception as e:
            logging.error(f"Error processing sample {true_folder}: {e}")

def main():
    base_path = "Favites_log"
    evaluator = StratusEvaluator()
    
    # Process all samples
    total_processed = 0
    for i in range(1, 501):
        true_folder = f"{base_path}/FAVITES_output_logSI_contemp_T200_N100_E1_{i}/"
        generated_folder = f"Inferred_tt/FAVITES_output_logSI_contemp_T200_N100_E1_{i}_tt.txt"
        
        if os.path.exists(generated_folder):
            evaluator.process_sample(true_folder, generated_folder)
            total_processed += 1
    
    # Save results
    evaluator.save_results()
    logging.info(f"Evaluation completed. Processed {total_processed} samples.")

if __name__ == "__main__":
    main()

