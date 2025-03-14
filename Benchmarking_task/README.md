# **Comparative Analysis of SMiTH Model on Favites Datasets**  

In this directory, we present a **comparative analysis** of the **SMiTH model** results on two simulated datasets: **Favites_log** and **Favites_exp**. The goal of this analysis is to evaluate the performance of **SMiTH** in reconstructing transmission networks compared to other models.

## **Simulated Datasets**
This section provides details about the simulated datasets used in our analysis. These datasets were generated using FAVITES (FrAmework for VIral Transmission and Evolution Simulation).

Intra-site evolution under the logistic coalescent model:  [Favites_log](https://uconn-my.sharepoint.com/:f:/g/personal/pavel_skums_uconn_edu/EkBG6Wu6y1dGpZFsFfWgJmUBaLBqtF6eRCBjlLW_BH-x1A?e=6wYFKX) 

Intra-site evolution under the exponential coalescent model: [Favites_exp](https://uconn-my.sharepoint.com/:f:/g/personal/pavel_skums_uconn_edu/EoTTY0RpWrZOgsCdK-qVHNgBx-DcakSYjwtzwVTicJ4Zmw?e=NHBQOd). 


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
