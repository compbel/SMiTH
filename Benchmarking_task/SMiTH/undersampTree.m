function [AMtree,patients,seqF] = undersampTree(filePhylo,fileSeq)

% This function was used for benchmarking of algorithms' performance under
% random sampling.
% Input:
% filePhylo: file with the phylogeny in newick format
% fileSeq: file with sequences in fasta format
% Output:
% AMtree: adjacency matrix of a phylogenetic tree for sequence subsample
% patients: list of site IDs of leafs of the phylogeny
% seqF: sequence subsample
% For more details on the subsampling method see the paper.
% To use the script, insert it to the script migratonSampler immediately
% before the line [AMtree, patients,~] = reduceTree(AMtree,patients);

phylotree = phytreeread(filePhylo);
[AMtree,patients] = phytree2graph(phylotree,delimeter,tokenPos);
tree = digraph(AMtree);   
seqF = fastaread(fileSeq);

names = get(phylotree,'LeafNames');
nSeq = length(names);
seqID = zeros(1,nSeq);
for i = 1:nSeq
    C = strsplit(names{i},'|');
    seqID(i) = str2num(C{1}(2:end));
end

seqIDfas = zeros(1,length(seqF));
for i = 1:length(seqF)
    C = strsplit(seqF(i).Header,'|');
    seqIDfas(i) = str2num(C{1}(2:end));
end
[~, indices] = ismember(seqID, seqIDfas);

u = unique(patients(patients~=0));
percSamp = rand(numel(u), 1); 

leafs = find(patients);          
[~, idx] = ismember(patients(leafs), u);

toRemove = leafs(rand(numel(leafs), 1) >= percSamp(idx));
tree = rmnode(tree, toRemove);
patients(toRemove) = []; 
seqF(indices(toRemove)) = [];
seqIDfas(indices(toRemove)) = [];

while any(outdegree(tree)==0 & patients'==0,1)
    leaf0 = find(outdegree(tree)==0 & patients'==0,1);
    tree = rmnode(tree, leaf0);
    patients(leaf0) = [];
end


while any(indegree(tree) == 1 & outdegree(tree) == 1)
    v = find(indegree(tree) == 1 & outdegree(tree) == 1, 1);
    p = predecessors(tree, v);  
    c = successors(tree, v); 
    tree = addedge(tree, p, c);
    tree = rmnode(tree, v);
    patients(v) = [];
end
badRoot = find(indegree(tree)==0 & outdegree(tree) == 1);
if ~isempty(badRoot)
    tree = rmnode(tree, v);
    patients(badRoot) = [];
end

AMtree = adjacency(tree);
