# **Comparative Analysis of SMiTH Model on Favites Datasets**  

In this directory, we present a **comparative analysis** of the **SMiTH model** results on two simulated datasets: **Favites_log** and **Favites_exp**. The goal of this analysis is to evaluate the performance of **SMiTH** in reconstructing transmission networks compared to other models.

## **Simulated Datasets**
This section provides details about the simulated datasets used in our analysis. These datasets were generated using FAVITES (FrAmework for VIral Transmission and Evolution Simulation).

Intra-site evolution under the logistic coalescent model:  [Favites_log](https://uconn-my.sharepoint.com/:u:/r/personal/marykafi_uconn_edu/Documents/FavitesDataset/Favites_log.zip?csf=1&web=1&e=x8HFsd) 

Intra-site evolution under the exponential coalescent model: [Favites_exp](https://uconn-my.sharepoint.com/:u:/r/personal/marykafi_uconn_edu/Documents/FavitesDataset/Favites_exp.zip?csf=1&web=1&e=PVR1wv). 

## **SMiTH scripts**
The folder **SMiTH** contains Matlab scripts used to run SMiTH on simulated data, calculate algorithm's accuracy and compare it with results of other tools. 

run_simulated_data.m: the main script that can be used to run SMiTH on all simulated datasets. Requires settings correct paths to input and output directories.

analyzePlotSimul.m: processing the algoritms' results and creating plots used in the paper

si_model_simul_rand_rates.m: implementation of a SI model with variable migration rates (used in benchmarking, the details can be found in the paper)

undersampTree.m: used to produce data under random site sampling.

## **Comparison Overview**
We benchamrked the **SMiTH tool** against the following transmission network inference tools:  
- **Cassiopeia**  
- **Machina**  
- **Phyloscanner**  
- **STraTUS**  
- **TNet**  


## **Directory Structure & Contents**
Each **model-specific directory** contains:  
**`/model_name/`** — A subdirectory for each model comparison, including:  
- **`run_model.sh` / `run_model.py`** — Code or scripts used to run the model on the dataset.  
- **Evaluation scripts** — Methods for analyzing the model’s performance.  



## **How to Use This Repository**
### **Navigate to the Model Directory**
Each model has its own folder where you can find scripts and evaluation files. 

### **Run the Model**
Follow the provided scripts and instructions to **replicate the execution** of the model on the dataset.
