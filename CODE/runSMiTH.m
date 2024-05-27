warning('off','all')
filePhylo = ['..' filesep 'input example' filesep 'small_example_input.csv'];
sampGenerator = @randTreePrefAttach;
nSamp = 100;
constr = 'compact'; 
timeLimit = NaN;
fileSeq = ['..' filesep 'input example' filesep 'small_example_sequence_data.fasta'];
% fileSeq = [];
delimeter = '|';
tokenPos = 2;


[migrSamp,objSamp,originSamp,consensus,siteList] = migrationSampler(filePhylo,sampGenerator,...
    nSamp,constr,timeLimit,fileSeq,delimeter,tokenPos);

% In these lines, we construct a maximal spanning tree of the consensus,
% and plot it
migrSpan = minspantree(graph(-consensus));
plot(migrSpan,'Layout','force','NodeLabel',string(siteList));
