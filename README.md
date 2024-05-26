# SMiTH
SMiTH (Sampling MIgration Trees using Homomorphisms) is a tool for inferring migration trees with expected properties.

## Algorithm
![alt text](/images/flow.png)

**A**: Input phylogenetic tree. 

**B**: Distribution of possible migration trees. Parallel rectangles depict a probability density function, with each rectangle's width proportional to the corresponding probability. In practice, the distribution is represented either by a random graph model or by a stochastic graph generation procedure. 

**C**: Candidate migration trees sampled from the distribution **B**. 

**D**: Homomorphisms from the phylogeny to three sampled trees. In the phylogenies, nodes are color-coded by their homomorphic images in a migration tree. The phylogeny layouts in the middle of each subfigure showcases how homomorphism transforms them into sampled trees.  

**E**: consensus solution derived from homomorphisms in **D**. The solution is shown as potential color distributions for the phylogeny's nodes (left) or as a graph where possible migration edges are weighted according to the number of supporting solutions (right), with the edge thickness indicating weight.


## Pre-requisites
   - Matlab
   - Gurobi

## Instructions

The project is run from the script ``runSMiTH.m``. In this script the user run the main function ``migrationSampler`` with necessary parameters.
``[migrSamp,objSamp,originSamp,consensus, siteList] = migrationSampler(filePhylo,sampGenerator,nSamp,constr,timeLimit,fileSeq,delimeter,tokenPos)``

### Input: required parameters
* ``filePhylo``:   file with the phylogenetic tree. Should consist of _N_ rows of the form:

   `p1` `s1`
   
   `p2` `s2`
   
   ...
   
   `pN` `sN`

   where
   - the _i_<sup>th</sup> row corresponds to the _i_<sup>th</sup> tree node;
   - `pi` is the parent of the node _i_ and ``si`` is an ID of the site corresponding to that node;
   - `pr=0` for the root node `r`, and `si=0` for internal nodes _i_.

* ``sampGeneratorf``   handle to the function sampling candidate migration trees from a particular distribution. Should have the tree size as a single argument. Current version of SMiTH package provides two predefined sampling function:
   - `randTreeUniform` for uniform sampling;
   - `randTreePrefAttach` for sampling via preferential attachment procedure.
     
  Users can provide their own custom sampling functions.
  
* ``nSamp``   number of candidate trees to be sampled
* ``constr``   structural constraints for inference of compatibility of a migration tree and a phylogeny.
     Possible values:
     - `'unconstrained'`;
     - `'convex'`;
     - `'convexMaxCompact'` (convex homomorphism with maximum compactness);
     - `'compact'`.

### Input: optional parameters
   * ``timeLimit``   time limit for running an ILP solver for the unconstrained homomorphism problem. Not relevant for other types of constraints, can be set as `[]`.
   * ``fileSeq``   fasta file with sequences. Should be specified only if genetic diversity of populations is used in calculations. It is assumed that for each sequence the ID of the population where it belongs is the part of its header. If you do not want to use diversity in the algorithm, set `divers = []`;
   * ``delimeter``   character used to specify the boundary between different tokens of a sequence header, with population ID being one of these tokens.
   * ``tokenPos``   the index of the token with the population ID.

### Output
 * ``migrSamp``   sample of migration trees compatible with the given phylogeny.
 * ``objSamp``   values of objective functions for sampled tree. 
 * ``originSamp``   origins (i.e., migration site corresponding to the root of the phylogeny) for sampled trees. Can be used to analyze migration directionality if needed.
 * ``consensus``   consensus matrix of the sample, i.e., `consensus(i,j)` is the frequency of sampled migration trees with vertices _i_ and _j_ being adjacent
 * ``siteList``   the list of migration sites in the same order as the migration tree vertices, i.e., `siteList(i)` is the site corresponding to the _i_<sup>th</sup> vertex of migration trees from `migrSamp`.

## Example: 
`filePhylo = ['input example' filesep 'input8.csv'];
sampGenerator = @randTreePrefAttach;
nSamp = 100;
constr = 'convexMaxCompact'; 
timeLimit = 600;
fileSeq = ['input example' filesep 'sequence_data8.fasta'];
delimeter = '|';
tokenPos = 2;

[migrSamp,objSamp,originSamp,consensus,siteList] = migrationSampler(filePhylo,sampGenerator,nSamp,constr,timeLimit,fileSeq,delimeter,tokenPos);`


## Citation:
TBA
