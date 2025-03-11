% runMigrationSampler.m - Read inputs from a config file

% Read input file (assumed format: key=value)
configFile = 'config.txt';
fileID = fopen(configFile, 'r');
configData = textscan(fileID, '%s %s', 'Delimiter', '=');
fclose(fileID);

% Convert inputs
inputMap = containers.Map(configData{1}, configData{2});

% Assign values (convert numbers where needed)
filePhylo = inputMap('filePhylo');
sampGenerator = @randTreePrefAttach;  % Assume this remains constant
nSamp = str2double(inputMap('nSamp'));
constr = inputMap('constr');
timeLimit = str2double(inputMap('timeLimit'));
perc = str2double(inputMap('perc'));
fileSeq = inputMap('fileSeq');
delimeter = inputMap('delimeter');
tokenPos = str2double(inputMap('tokenPos'));

% Run main function
[originSamp,consensus,siteList] = migrationSamplerReduced(...
    filePhylo,sampGenerator,nSamp,constr,timeLimit,perc,fileSeq,delimeter,tokenPos);

% Display results
disp('Migration Sampling Completed!');
disp(migrSamp);