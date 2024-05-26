warning('off','all')
filePhylo = ['input example' filesep 'input8.csv'];
sampGenerator = @randTreePrefAttach;
nSamp = 100;
constr = 'convexMaxCompact'; 
timeLimit = 600;
% fileSeq = ['input example' filesep 'sequence_data8.fasta'];
fileSeq = [];
delimeter = '|';
tokenPos = 2;


[migrSamp,objSamp,originSamp,consensus, siteList] = migrationSampler(filePhylo,sampGenerator,...
    nSamp,constr,timeLimit,fileSeq,delimeter,tokenPos);

% In these lines, we construct a maximal spanning tree of the consensus,
% and plot it
migrSpan = minspantree(graph(-consensus));
plot(migrSpan,'Layout','force','NodeLabel',string(siteList));
