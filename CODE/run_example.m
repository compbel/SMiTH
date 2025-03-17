clear;
% filePhylo = ['input example' filesep 'tree16.csv'];
filePhylo = ['input example' filesep 'tree16.tre'];
sampGenerator = @randTreePrefAttach;
nSamp = 10000;
constr = 'convexMaxCompact'; 
fileSeq = ['input example' filesep 'sequence_data16.fasta'];
delimeter = '|';
tokenPos = 2;
timeLimit = 600;

[migrSamp,objSamp,originSamp,consensus, siteList] = migrationSampler(filePhylo,sampGenerator,...
    nSamp,constr,timeLimit,fileSeq,delimeter,tokenPos);

% In these lines, we constract a maximal spanning tree of the consensus,
% and plot it
migrSpan = minspantree(graph(-consensus));
plot(migrSpan,'Layout','force','NodeLabel',string(siteList));