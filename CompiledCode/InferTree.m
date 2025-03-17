% Read input file
%configFile = 'config.txt';
function [originSamp,consensus,siteList] = InferTree(configFile)
fileID = fopen(configFile, 'r');
configData = textscan(fileID, '%s %s', 'Delimiter', '=');
fclose(fileID);

% Convert inputs
inputMap = containers.Map(configData{1}, configData{2});

% Assign values (convert numbers where needed)
dir_name = inputMap('dir_name');
filePhylo = inputMap('filePhylo');
sampGenerator = @randTreePrefAttach;  % Assume this remains constant
nSamp = str2double(inputMap('nSamp'));
constr = inputMap('constr');
timeLimit = str2double(inputMap('timeLimit'));
perc = [];
fileSeq = inputMap('fileSeq');
delimeter = inputMap('delimeter');
tokenPos = str2double(inputMap('tokenPos'));






% main function generating a list of trees compatible with the given
% phylogeny and sampled from the given random tree distribution
% Input:
% Required parameters:
% (*) filePhylo - file with the phylogenetic tree. Should consist of N rows
% of the form
% p1 s1
% p2 s2
% ...
% pN sN
% where the ith row corresponds to the ith tree node, pi is the
% parent of the node i and pi is an id of the site corresponding to that
% node. pr=0 for the root node r, and si=0 for internal nodes i.
% (*) sampGenerator - handle to the function sampling candidate migration trees from a
% particular distribution. Should have the tree size as a single argument.
% Current version of SMiTH package provides two predefined sampling
% function: randTreeUniform for uniform sampling, and randTreePrefAttach
% for sampling via preferential attachment procedure. Users can proved
% their own custom sampling functions.
% (*) nSamp - number of candidate trees to be sampled.
% (*) constr - structural constraints for inference of compatibility of a
% migration tree and a phylogeny. Possible values: 'unconstrained', 
% 'convex', 'convexMaxCompact' (convex homomorphism with maximum
% compactness), 'compact'.
% Optional parameters:
% (*) timeLimit: time limit for running an ILP solver for the unconstrained
% homomorphism problem. Not relevant for other types of constraints, can be
% set as [];
% (*) perc: percentile of sampled migration trees ranked by their objective values
% that should be used in an output. Can be set as []; in that case the
% entire sample is used;
% (*) fileSeq: fasta file with sequences. Should be specified only if  genetic 
% diversity of populations is used in calculations. It is assumed that for each
% sequence the id of the population where it belongs is the part of its header. 
% If you do not want to use diversity in the algorithm, set divers = [];
% (*) delimeter: character used to specify the boundary between different tokens 
% of a sequence header, with population id being one of these tokens.
% (*) tokenPos: the index of the token with the population id.
% Output:
% (*) originSamp - frequency of origins (i.e. migration site corresponding to the root of the phylogeny)
% for sampled trees. Can be used to analyze migration directionality if
% needed.
% (*) consensus - consensus matrix of the sample, i.e. consensus(i,j) is
% the frequency of sampled migration trees with vertices i and j being
% adjacent
% (*) siteList - the list of migration sites in the same order as the
% migration tree vertices, i.e. siteList(i) is the site corresponding to
% the ith vertex of migration trees from migrSamp.

contrUnconstr = 0;
treeData = readmatrix(filePhylo);
AMtree = zeros(size(treeData,1),size(treeData,1));
patients = zeros(1,2*size(treeData,1)-1);
for i = 1:size(treeData,1)
    if treeData(i,1)~=0
        AMtree(treeData(i,1),i) = 1;
    end
    patients(i) = treeData(i,2);
end

[AMtree, patients,~] = reduceTree(AMtree,patients);
patientList = sort(unique(patients));
patientList = patientList(2:end);
nPat = length(patientList);

if ~isempty(fileSeq)
    divers = getDiversity(fileSeq,delimeter,tokenPos);
else
    divers = [];
end

[traits,~] = pat2traits(patients,patientList);

if strcmp(constr,'convex') | strcmp(constr,'compact') | strcmp(constr,'convexMaxCompact') | contrUnconstr
    [AMtree,traits] = contractMultLabel(AMtree,traits,nPat);
    [AMtree,traits] = renameNodeTraits(AMtree,traits,nPat);
end


tree = digraph(AMtree>0);
AMSamp = cell(1,nSamp);
objSamp = -ones(1,nSamp);
originSamp = -ones(1,nSamp);

parfor s = 1:nSamp
    disp(['Sampled tree ' num2str(s)]);
    cand = sampGenerator(nPat);
    if strcmp(constr,'unconstrained')
        [ishom,AMTNinferCurr,obj,origin] = checkHomomorUncons(tree,cand,traits,divers,timeLimit);
    end
    if strcmp(constr,'convex')
        toOptimizeComp = 0;
        [ishom,AMTNinferCurr,obj,origin] = checkHomomorConvex(tree,cand,traits,toOptimizeComp,divers);
    end
    if strcmp(constr,'convexMaxCompact')
        toOptimizeComp = 1;
        [ishom,AMTNinferCurr,obj,origin] = checkHomomorConvex(tree,cand,traits,toOptimizeComp,divers);
    end
    if strcmp(constr,'compact')
        [ishom,AMTNinferCurr,obj,origin] = checkHomomorCompact(tree,cand,traits);
    end
    if ishom
        objSamp(s) = obj;
        originSamp(s) = origin;
        AMsamp{s} = AMTNinferCurr;
    end
end

ind = (objSamp ~= -1);
AMsamp = AMsamp(ind);
objSamp = objSamp(ind);
originSamp = originSamp(ind);

if ~isempty(perc)
    ind = (objSamp >= perc);
    AMsamp = AMsamp(ind);
    originSamp = originSamp(ind);
end


histEdges = 0.5:1:(nPat+0.5);
originSamp = histcounts(originSamp, histEdges);
consensus = sum(cat(3, AMsamp{:}), 3)/length(AMsamp);
siteList = patientList;
% Display results


% Create output directory
outputDir = fullfile('output', dir_name);
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% Save files
writematrix(originSamp, fullfile(outputDir, 'originSamp.csv'));
writematrix(consensus, fullfile(outputDir, 'consensus.csv'));
writematrix(siteList, fullfile(outputDir, 'siteList.csv'));

disp(['Data saved in: ', outputDir]);
disp('Migration Sampling Completed!');
end