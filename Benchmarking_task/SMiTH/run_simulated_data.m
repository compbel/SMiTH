clear;
warning('off')
prior = 'SF'; %Prufer %True
coal = 'exp';
sampGenerator = @randTreePrefAttach;
delimeter = '|';
tokenPos = 2;
timeLimit = 600;
nSamp = 150000;
constr = 'convexMaxCompact'; 
nTests = 500;


sens = zeros(nTests,1);
spec = zeros(nTests,1);
fscore = zeros(nTests,1);
runtime = zeros(nTests,1);

resFolder = ['results' filesep prior '_' coal];
if ~exist(resFolder, 'dir')
    mkdir(resFolder);
end
dataFolder = 'Favites example';

for t = 1:nTests
    resFile = [resFolder filesep 'res_t' num2str(t) '.mat'];
    favOutFolder = [dataFolder filesep 'FAVITES_output_' coal 'SI_contemp_T200_N100_E1_' int2str(t)];

    fileCN = [favOutFolder filesep 'contact_network.txt'];
    fileSeq = [favOutFolder filesep  'error_free_files' filesep 'sequence_data.fasta'];
    fileTN = [favOutFolder filesep 'error_free_files' filesep 'transmission_network.txt'];
    filePhylo = [favOutFolder filesep 'error_free_files' filesep 'phylogenetic_trees' filesep 'tree_0.time.tre'];
    
    dataCN = readlines(fileCN);
    C = strsplit(dataCN(1));
    nCN = 1;
    while strcmp(C(1),'NODE')
        nCN = nCN + 1;
        C = strsplit(dataCN(nCN));
    end
    nCN = nCN - 1;

    AMTNtrue = zeros(nCN,nCN);
    dataTN = readlines(fileTN);
    for i = 1:size(dataTN,1)-1
        C = strsplit(dataTN(i));
        if strcmp(C(1),'None')
            source = str2num(C(2))+1;
        else
            u = str2num(C(1))+1;
            v = str2num(C(2))+1;
            AMTNtrue(u,v) = 1;
        end
    end    
    AMTNtrue = AMTNtrue - diag(diag(AMTNtrue));

    nVertTest = sum(sum(AMTNtrue))+1;
    if ~ismember(nVertTest,5:30)
        continue;
    end
    
    tic
    [migrSamp,objSamp,originSamp,consensus, siteList] = migrationSampler(filePhylo,sampGenerator,...
        nSamp,constr,timeLimit,fileSeq,delimeter,tokenPos);
    timeSamp = toc;
    AMTNtrue =  AMTNtrue(siteList,siteList);
    save(resFile,'migrSamp','AMTNtrue','objSamp','timeSamp');
end