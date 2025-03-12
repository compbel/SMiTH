# SMiTH
SMiTH (**S**ampling **Mi**gration **T**rees using **H**omomorphisms) is a tool for inferring migration trees with expected properties.

## Algorithm
![alt text](/images/flow.png)

**A**: Input phylogenetic tree. 

**B**: Distribution of possible migration trees. Parallel rectangles depict a probability density function, with each rectangle's width proportional to the corresponding probability. In practice, the distribution is represented either by a random graph model or by a stochastic graph generation procedure. 

**C**: Candidate migration trees sampled from the distribution **B**. 

**D**: Homomorphisms from the phylogeny to three sampled trees. In the phylogenies, nodes are color-coded by their homomorphic images in a migration tree. The phylogeny layouts in the middle of each subfigure showcases how homomorphism transforms them into sampled trees.  

**E**: consensus solution derived from homomorphisms in **D**. The solution is shown as potential color distributions for the phylogeny's nodes (left) or as a graph where possible migration edges are weighted according to the number of supporting solutions (right), with the edge thickness indicating weight.


## Running SMiTH Without a MATLAB License  
This guide provides step-by-step instructions for running the compiled MATLAB application inferTree on a machine without MATLAB installed. The application was created using MATLAB Compiler, and it requires the MATLAB Runtime to execute.
Additionally, Gurobi is needed only when using the 'unconstrained' mode for inferring the compatibility of a migration tree and a phylogeny.

### **Prerequisites**  
#### **1. MATLAB Runtime**
- The target machine **must have MATLAB Runtime installed**.  
- MATLAB Runtime is **free** and does **not require a MATLAB license**.  

#### **2. Operating System Compatibility**  
- Ensure the **target machine** runs the same operating system (Windows, macOS, or Linux) that the application was compiled for.  


### **Step 1: Download and Install MATLAB Runtime**  
MATLAB Runtime is required to run the compiled application.  

1. **Go to the MATLAB Runtime download page**:  
   - [MATLAB Runtime Download Page](https://www.mathworks.com/products/compiler/matlab-runtime.html)  
2. **Download the correct MATLAB Runtime version** that matches the version used to compile the application.  
3. **Install MATLAB Runtime**:
   - **Windows**: Run the installer executable (`.exe`) and follow the prompts.  
   - **Linux/macOS**: Open a terminal and run:  
     ```sh
     ./install -mode silent -agreeToLicense yes
     ```

### **Step 2: Download the Compiled Application**  
Ensure you have all the necessary compiled application files. These should include:  

#### **Executable File**  
- **Windows:** `CODE/runInferTree.exe`  
- **macOS:** `CODE/runInferTree.app.zip` *(Extract before running)*  
- **Linux:** `CODE/runInferTree.sh`  

#### **Shell Script to Launch the Application**  
- `CODE/run_runInferTree.sh` *(Required for macOS/Linux users)*  

#### **Additional Input Files & Dependencies**  
- `CODE/config.txt` *(Configuration file for setting input parameters)*  
- Any required **data files** for the program  




### **Step 3: Configure Input Parameters (`config.txt`)**  
Before running the program, update `CODE/config.txt` to define the required parameters for runMigrationSamplerReduced. This file specifies key input settings, such as the phylogenetic tree file, sampling method, and structural constraints. Properly configuring this file ensures that the program runs with the correct input data.
 
### **Required Parameters**
* ``filePhylo``:  csv-file with the phylogenetic tree. It must consist of _N_ rows and 2 columns `parent ID` and `label` (color of this node / ID of the site corresponding to this node):
   - the _i_<sup>th</sup> row corresponds to the _i_<sup>th</sup> tree node;
   - `pi` is the parent of the node _i_ and ``si`` is an ID of the site corresponding to that node;
   - `pr=0` for the root node `r`, and `si=0` for internal nodes _i_.

   `p1` `s1`
   
   `p2` `s2`
   
   `...`
   
   `pN` `sN`
  
     See some examples in the folder `input example`.

* ``sampGenerator``   handle to the function sampling candidate migration trees from a particular distribution. Should have the tree size as a single argument. Current version of SMiTH package provides two predefined sampling function:
   - `@randTreeUniform` for uniform sampling;
   - `@randTreePrefAttach` for sampling via preferential attachment procedure.
     
  Users can provide their own custom sampling functions.
  
* ``nSamp``   number of candidate trees to be sampled.
* ``constr``   structural constraints for inference of compatibility of a migration tree and a phylogeny.
     Possible values:
     - `'unconstrained'`
     - `'convex'`
     - `'convexMaxCompact'` (convex homomorphism with maximum compactness)
     - `'compact'`

### Input: optional parameters
   * ``timeLimit``   time limit for running an ILP solver for the unconstrained homomorphism problem. If ``constr='unconstrained'``, set it to number of second Gurobi solves ILP; if unlimited, set `timeLimit = NaN`. If `constr` is not `'unconstrained'`, set `timeLimit = NaN`.

   * ``perc``   percentile of sampled migration trees ranked by their objective values that should be used in an output. Can be set as ``perc = []``. in that case the entire sample is used;

   * ``fileSeq``   fasta file with sequences. Should be specified only if genetic diversity of populations is used in calculations. It is assumed that for each sequence the ID of the population where it belongs is the part of its header. If you do not want to use diversity in the algorithm, set `fileSeq = []`.
   * ``delimeter``   character used to specify the boundary between different tokens of a sequence header in the fasta file, with population ID being one of these tokens. It is needed only when `fileSeq` is specified. E.g., for the sequence header `>N614|56|100.0`, it is a vertical bar `|`. If you do not use diversity in the algorithm, set `delimeter = NaN`.
   * ``tokenPos``   the index of the token of a sequence header in the fasta file with the population ID.  E.g., for the sequence header `>N614|56|100.0`, it is equal `2` (we want to separate `56`). If you do not use diversity in the algorithm, set `tokenPos = NaN`.
### Example
```
filePhylo=myPhyloFile.txt
nSamp=200
constr=convexMaxCompact
timeLimit=500
perc=100
fileSeq=sequence_data8.fasta
delimeter=|
tokenPos=2
```

### **Expected Output**
This part describes the **output files/variables** that the user will obtain **after running `runMigrationSamplerReduced`**.

 * ``originSamp``   origins (i.e., migration site corresponding to the root of the phylogeny) for sampled trees. Can be used to analyze migration directionality if needed.
 * ``consensus``   consensus matrix of the sample, i.e., `consensus(i,j)` is the frequency of sampled migration trees with vertices _i_ and _j_ being adjacent.
 * ``siteList``   the list of migration sites in the same order as the migration tree vertices, i.e., `siteList(i)` is the site corresponding to the _i_<sup>th</sup> vertex of migration trees from `migrSamp`.




### **Step 4: Run the Application**  
#### **Windows**  
1. Open **Command Prompt** (`cmd`).  
2. Navigate to the directory containing the compiled files:  
```sh
cd CODE
runInferTree.exe
```

#### **Linux/MacOS**  
1. Open **Command Prompt** (`cmd`).  
2. Navigate to the directory containing the compiled files:  
```sh
cd Compiler
chmod +x run_runInferTree.sh
./run_runInferTree.sh
```


## Running SMiTH Using MATLAB License 
This guide provides step-by-step instructions for running the SMiTH code with a MATLAB license. The program requires MATLAB to be installed, and Gurobi is needed only when using the 'unconstrained' mode for inferring the compatibility of a migration tree and a phylogeny.

The project is executed using the runSMiTH.m script, where the user should call the main function migrationSampler with the required parameters.

```[migrSamp,objSamp,originSamp,consensus,siteList] = migrationSampler(filePhylo,sampGenerator,nSamp,constr,timeLimit,fileSeq,delimeter,tokenPos)```

### Input: required parameters
* ``filePhylo``:  csv-file with the phylogenetic tree. It must consist of _N_ rows and 2 columns `parent ID` and `label` (color of this node / ID of the site corresponding to this node):
   - the _i_<sup>th</sup> row corresponds to the _i_<sup>th</sup> tree node;
   - `pi` is the parent of the node _i_ and ``si`` is an ID of the site corresponding to that node;
   - `pr=0` for the root node `r`, and `si=0` for internal nodes _i_.

   `p1` `s1`
   
   `p2` `s2`
   
   `...`
   
   `pN` `sN`
  
     See some examples in the folder `input example`.

* ``sampGenerator``   handle to the function sampling candidate migration trees from a particular distribution. Should have the tree size as a single argument. Current version of SMiTH package provides two predefined sampling function:
   - `@randTreeUniform` for uniform sampling;
   - `@randTreePrefAttach` for sampling via preferential attachment procedure.
     
  Users can provide their own custom sampling functions.
  
* ``nSamp``   number of candidate trees to be sampled.
* ``constr``   structural constraints for inference of compatibility of a migration tree and a phylogeny.
     Possible values:
     - `'unconstrained'`
     - `'convex'`
     - `'convexMaxCompact'` (convex homomorphism with maximum compactness)
     - `'compact'`

### Input: optional parameters
   * ``timeLimit``   time limit for running an ILP solver for the unconstrained homomorphism problem. If ``constr='unconstrained'``, set it to number of second Gurobi solves ILP; if unlimited, set `timeLimit = NaN`. If `constr` is not `'unconstrained'`, set `timeLimit = NaN`.
   * ``fileSeq``   fasta file with sequences. Should be specified only if genetic diversity of populations is used in calculations. It is assumed that for each sequence the ID of the population where it belongs is the part of its header. If you do not want to use diversity in the algorithm, set `fileSeq = []`.
   * ``delimeter``   character used to specify the boundary between different tokens of a sequence header in the fasta file, with population ID being one of these tokens. It is needed only when `fileSeq` is specified. E.g., for the sequence header `>N614|56|100.0`, it is a vertical bar `|`. If you do not use diversity in the algorithm, set `delimeter = NaN`.
   * ``tokenPos``   the index of the token of a sequence header in the fasta file with the population ID.  E.g., for the sequence header `>N614|56|100.0`, it is equal `2` (we want to separate `56`). If you do not use diversity in the algorithm, set `tokenPos = NaN`.

### Output
 * ``migrSamp``   sample of migration trees compatible with the given phylogeny.
 * ``objSamp``   values of objective functions for sampled tree. 
 * ``originSamp``   origins (i.e., migration site corresponding to the root of the phylogeny) for sampled trees. Can be used to analyze migration directionality if needed.
 * ``consensus``   consensus matrix of the sample, i.e., `consensus(i,j)` is the frequency of sampled migration trees with vertices _i_ and _j_ being adjacent.
 * ``siteList``   the list of migration sites in the same order as the migration tree vertices, i.e., `siteList(i)` is the site corresponding to the _i_<sup>th</sup> vertex of migration trees from `migrSamp`.

### Example
```filePhylo = ['input example' filesep 'input8.csv'];
sampGenerator = @randTreePrefAttach;
nSamp = 100;
constr = 'convexMaxCompact'; 
timeLimit = 600;
fileSeq = ['input example' filesep 'sequence_data8.fasta'];
delimeter = '|';
tokenPos = 2;

[migrSamp,objSamp,originSamp,consensus,siteList] = migrationSampler(filePhylo,sampGenerator,...
    nSamp,constr,timeLimit,fileSeq,delimeter,tokenPos);
```
## **Bechmarking Study**
The Bechmarking_task directory provides a comparative analysis of the SMiTH model performance on two simulated datasets, Favites_log and Favites_exp, assessing its ability to reconstruct transmission networks compared to other models. Both datasets simulate disease transmission and viral evolution using a Barabási-Albert contact network with 100 nodes and the Generalized Epidemic Modeling Framework (GEMF). The datasets differ in epidemic growth models, with Favites_log following a logistic growth pattern (gradual spread and plateau) and Favites_exp using an exponential growth model (rapid outbreak expansion).

The analysis compares SMiTH to various inference models, including Cassiopeia, Machina, Phyloscanner, STraTUS, TNet, and TransPhyloMulti. Each model-specific directory contains scripts for running the model, evaluating its performance, and comparing results against the ground truth transmission network using specificity, sensitivity, and F1-score metrics. The repository is structured to allow easy navigation, execution, and evaluation of different models on the Favites datasets.

## Citation
TBA
