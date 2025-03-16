clear;


fold = 'secondary data';
dataset = 'mcpherson_2016';
prior = 'uniform';
patIDs = [1 4 7];
nSamp = 50000;
objDeg = 0;
divers = [];

%% 

for i = patIDs
    outFile = ['results' filesep 'realCancer' filesep 'res_' dataset '_' int2str(i) '.mat']; %output file for sampling results
    dataFile = [fold filesep dataset '_nodesdata_' int2str(i) '.csv'];
    data = readtable(dataFile);
    sites = unique(data.anatomicalSite);
    sites(strcmp(sites,'NA')) = [];
    nPat = length(sites);
    
        if i == 7
            treeResolved = [16 1 2; 17 16 3; 18 17 4; 19 18 5; 20 19 6; 21 20 7; 22 21 8; 23 22 9; 24 23 10;...
            25 24 11; 26 14 15; 27 26 13; 28 27 12; 29 25 28];
            % Polytomies are resolved in such a way so that minimal migration/co-migration numbers
            % reported by Machina in supplement achievable. Everything else is random
            organs = {'LOv','RUt','RPv','Bwl','Bwl','RUt','LOv','ROv','ROv','LOv','Bwl','RUt','RUt','Brn','Bm'};
            toPermut = [2 4; 5 8; 9 11];
        end
        if i == 1
            treeResolved = [20 1 2; 21 20 3; 22 21 4; 23 22 5; 24 23 29; 25 6 7; 26 25 8; 27 26 9; 28 27 10; 29 28 11;...
                30 12 13; 31 30 14; 32 15 16; 33 31 32; 34 33 17; 35 34 18; 36 35 19; 37 24 36];
            organs = {'Om','SBwl','LFTB','LOv','ApC','SBwl','SBwl','LOv','LFTB','ApC','Om','ROv','ROv','LOv','ROv','ROv','ROv','LOv','RFTA'};
            toPermut = [1 5; 8 11];
        end
        if i == 3
            treeResolved = [61 1 60; 60 2 59; 59 3 58; 58 4 57; 57 5 56; 56 6 55; 55 7 54; 54 8 53; 51 9 10;...
                52 51 11; 53 52 50; 50 31 49; 49 30 48; 48 29 47; 47 28 46; 46 27 45; 45 26 44; 44 25 43; ...
                43 24 42; 42 41 38; 39 12 13; 40 14 39; 41 15 40; 38 23 37; 37 22 36; 36 21 35; 35 20 34; ...
                34 19 33; 33 18 32; 32 16 17];
            
            organs = {'LOv','CDSB','ClnE','Adnx','LFTC','Om','ROv','RFTA','ClnE','CDSB','Om','LFTC','Adnx','Om','LOv',...
                'ROv','ROv','CDSB','RFTA','ClnE','LFTC','Om','LOv','LOv','ROv','RFTA','LFTC','Adnx','CDSB','ClnE','Om'};
            toPermut = [1 8; 9 11; 12 15; 16 23; 24 31];
        end
        n = 2*size(treeResolved,1)+1;
        AM = zeros(n,n);
        for j = 1:size(treeResolved,1)
            AM(treeResolved(j,1),treeResolved(j,2))=1;
            AM(treeResolved(j,1),treeResolved(j,3))=1;
        end        
        traits = zeros(1,length(AM));
        for j = 1:length(sites)
            ind = strcmp(organs,sites{j});
            traits(ind) = j;
        end

    tree = digraph(AM);
    node_labels = cell(1,numnodes(tree));
    for j = 1:length(node_labels)
        if traits(j) > 0
            node_labels{j} = sites{traits(j)};
        else
            node_labels{j} = '';
        end
    end
    figure
    p = plot(tree,'Layout','layered','NodeLabel',node_labels,'LineWidth',2);
    p.MarkerSize = 10; 
    p.NodeFontSize = 14;
    p.NodeFontWeight = 'bold';
    treeUndir = graph(AM+AM' > 0);

    predec = findPredecessors(tree,nPat);
    AMTNinferSamp = cell(1,nSamp);
    AMTNinferSampDir = cell(1,nSamp);
    objLPSamp = -ones(1,nSamp);
    objParsSamp = -ones(1,nSamp);
    originSamp = zeros(1,nSamp);
    parfor s = 1:nSamp
        s
        cand = randTreeUniform(nPat);
        traitsCurr = traits;
        for j = 1:size(toPermut,1)
            tr = traitsCurr(toPermut(j,1):toPermut(j,2));
            tr = tr(randperm(length(tr)));
            traitsCurr(toPermut(j,1):toPermut(j,2)) = tr;
        end
        try
            [ishom,map,invmap,AMTNinferCurr,objLP,objPars,suff,AMTNinferCurrDir,origin] = checkHomomorUnconstrMult(treeUndir,tree,cand,traitsCurr,objDeg,divers);
        catch
            continue;
        end
        AMTNinferSamp{s} = AMTNinferCurr;
        AMTNinferSampDir{s} = AMTNinferCurrDir;
        objLPSamp(s) = objLP;
        objParsSamp(s) = objPars;
        originSamp(s) = origin;
    end
    save(outFile,'AMTNinferSamp','objLPSamp','objParsSamp','AMTNinferSampDir','originSamp');
end

%% 
for i = patIDs
    outFile = ['results' filesep 'realCancer' filesep 'res_' dataset '_' int2str(i) '.mat'];
    load(outFile);
    dataFile = [fold filesep dataset '_nodesdata_' int2str(i) '.csv'];
    data = readtable(dataFile);    sites = unique(data.anatomicalSite);
    sites(strcmp(sites,'NA')) = [];

    ind = cellfun(@isempty,AMTNinferSamp);
    AMTNinferSamp = AMTNinferSamp(~ind);
    AMTNinferSampDir = AMTNinferSampDir(~ind);
    objParsSamp = objParsSamp(~ind);
    objLPSamp = objLPSamp(~ind);
    originSamp = originSamp(~ind);

    
    [originCounts, ~] = histcounts(originSamp, 'BinMethod', 'integers');
    originCounts = originCounts/sum(originCounts);

    figure
    boxchart(originSamp,objParsSamp,'GroupByColor',originSamp,...
        "BoxWidth",6,'BoxFaceAlpha',0.5,'BoxMedianLineColor','k')
    title('Patient 3','FontSize', 18);
    lgd = legend(sites,'Location','southoutside','FontSize',11,'FontWeight','bold');
    lgd.NumColumns = length(sites);
    ylabel('Migration number','FontSize', 14,'FontWeight','bold');
    ax = gca;
    ax.FontWeight = 'bold';
    set(gca,'XTickLabel',[]);
    grid on
    exportgraphics(gcf,['figures' filesep 'migNum_distrib' int2str(i) '.png'],'Resolution',600)


    kruskalwallis(objParsSamp,originSamp);
    allKS = zeros(length(sites),length(sites));
    allMW = zeros(length(sites),length(sites));
    for u = 1:length(sites)
        for v = 1:length(sites)
            p = ranksum(objParsSamp(originSamp == u),objParsSamp(originSamp == v));
            allMW(u,v) = p;
        end
    end
end

%% 
function [ishom,map,invmap,AM,objLP,objPars,suff,AMdir,origin] = checkHomomorUnconstrMult(treeUndir,tree,cand,traits,objDeg,divers)

timeLimit = 600;
eps = 0.001;
n = numnodes(cand);

indTraits = (traits > 0);
degCand = degree(cand);
degTree = divers;

DtreeAll = distances(treeUndir);
Dtree = zeros(n,n);
for i = 1:n
    for j = 1:n
        ind1 = (traits == i);
        ind2 = (traits == j);
        Dtree(i,j) = min(min(DtreeAll(ind1,ind2)));
    end
end
Dcand = distances(cand);

ishom = false;
map = [];
invmap = [];
objLP = -Inf;
objPars = Inf;
suff = -1;

id = randi(1000000);
filename = ['ILP_' int2str(id) '.lp'];
fileID = fopen(filename,'w');


    fprintf(fileID,'%s\n','maximize');    
    for i = 1:n
        for j = 1:n
            Xij = ['X' int2str(i) ',' int2str(j)];
            % fprintf(fileID,'%s',['+ ' num2str(abs(eccTree(i)-eccCand(j))) ' ' Xij ' ']);  
            if objDeg == 1
                fprintf(fileID,'%s',['+ ' num2str(degTree(i)*degCand(j)) ' ' Xij ' ']);  
            else
                fprintf(fileID,'%s',['+ ' Xij ' ']);  
            end
        end
    end




fprintf(fileID,'\n%s\n',' st');

for i = 1:n
    for j = 1:n
        Xij = ['X' int2str(i) ',' int2str(j)];
        fprintf(fileID,'%s',['+ ' Xij ' ']);  
    end
    fprintf(fileID,'%s\n',['= 1']);
end

for j = 1:n
    for i = 1:n
        Xij = ['X' int2str(i) ',' int2str(j)];
        fprintf(fileID,'%s',['+ ' Xij ' ']);  
    end
    fprintf(fileID,'%s\n',['= 1']);
end

for i = 1:n
    for j = (i+1):n
        for u = 1:n
            for v = 1:n
                if Dtree(i,j) < Dcand(u,v)
                    Xiu = ['X' int2str(i) ',' int2str(u)];
                    Xjv = ['X' int2str(j) ',' int2str(v)];
                    fprintf(fileID,'%s\n',[Xiu ' + ' Xjv ' <= 1']);
                end
            end
        end
    end
end

for i = 1:n
    for j = (i+1):n
        for a = 1:size(cand.Edges.EndNodes,1)
            u = cand.Edges.EndNodes(a,1);
            v = cand.Edges.EndNodes(a,2);
            Xiu = ['X' int2str(i) ',' int2str(u)];
            Xiv = ['X' int2str(i) ',' int2str(v)];
            Xju = ['X' int2str(j) ',' int2str(u)];
            Xjv = ['X' int2str(j) ',' int2str(v)];
            Yij = ['Y' int2str(i) ',' int2str(j)];
            fprintf(fileID,'%s\n',[Xiu ' + ' Xiv ' + ' Xju ' + ' Xjv ' - ' Yij ' <= 1']);   
        end
    end
end

for i = 1:n
    for j = (i+1):n
        Yij = ['Y' int2str(i) ',' int2str(j)];
        fprintf(fileID,'%s',['+ ' Yij ' ']);  
    end
end
fprintf(fileID,'%s\n',['= ' num2str(n-1)]);


fprintf(fileID,'%s\n','binaries');
for i = 1:n
    for j = 1:n
        Xij = ['X' int2str(i) ',' int2str(j)];
        fprintf(fileID,'%s\n',Xij);
    end
end


for i = 1:n
    for j = (i+1):n
        Yij = ['Y' int2str(i) ',' int2str(j)];
        fprintf(fileID,'%s\n',Yij);
    end
end


fprintf(fileID,'%s\n','end');
fclose(fileID);

clear model;
model = gurobi_read(filename);

params.outputflag = 1;
params.timelimit = timeLimit;
result = gurobi(model, params);
delete(filename);

AM = [];
AMdir = [];
origin = -Inf;
if isfield(result,'x')  
    ind = contains(model.varnames,'X');
    X = result.x(ind);
    M = (reshape(X,n,n))';
    M = round(M);
    ishom = 1;
    map = zeros(1,n);
    invmap = zeros(1,n);
    for i = 1:n
        map(i) = find(M(i,:) >= 1-eps);
        invmap(i) = find(M(:,i));
    end
    AM = zeros(n,n);
    for i = 1:n
        for j = (i+1):n
            Yij = ['Y' int2str(i) ',' int2str(j)];
            ind = find(strcmp(Yij,model.varnames));
            AM(i,j) = result.x(ind);
            AM(j,i) = result.x(ind);
        end
    end
    objLP = result.objval;
else
    return;
end

root = find(indegree(tree) == 0);
leafs = (outdegree(tree) == 0);
dfsorder = flip(dfsearch(tree,root))';
AMrefl = AM + eye(n,n);


possLabels = cell(1,numnodes(tree));
for v = dfsorder
    if leafs(v)
        possLabels{v} = traits(v);
    else
        child = successors(tree,v);
        ind = true(n,1);
        for i = 1:length(child)
            ind = ind & (sum(AMrefl(:,possLabels{child(i)}),2)>0); 
        end
        possLabels{v} = find(ind);

        if isempty(possLabels{v})
            ishom = false;
            map = [];
            invmap = [];
            objLP = -Inf;
            objPars = Inf;
            suff = 0;
            AM = [];
            return;
        end
    end
end
suff = 1;
labels = zeros(1,numnodes(tree));
dfsorder = flip(dfsorder);
objPars = 0;
AMdir = zeros(n,n);
for v = dfsorder
    if v == root
        if ~isempty(divers)
            i = find(divers(possLabels{v})==max(divers(possLabels{v})),1);
            labels(v) = possLabels{v}(i);
        else
            labels(v) = possLabels{v}(randi(numel(possLabels{v})));
        end
        origin = labels(v);
        continue;
    end
    p = predecessors(tree,v);
    ind = find(AMrefl(labels(p),possLabels{v}));
    lb = possLabels{v}(ind);
    if ismember(labels(p),lb)
        labels(v) = labels(p);
    else
        if ~isempty(divers)
            i = find(divers(lb) == max(divers(lb)),1);
            labels(v) = lb(i);
        else
            labels(v) = lb(randi(numel(lb)));
        end
        objPars = objPars + 1;
        AMdir(labels(p),labels(v)) = 1;
    end
end
end

