clear;

fold = 'secondary data';
outbreaks = {'AI','AW'};
sources = [4 2];
nOut = length(outbreaks);
nSamp = 100000;
toVis = 1;
objDeg = 1;
nodeLabels = 'traits';


%% 
for i = 1:nOut
    % i = 3
    outb_name = outbreaks{i};
    outFile = ['results' filesep 'realHCV' filesep 'res_' outb_name '.mat'];
    filePhylo = [fold filesep 'RAxML_bestTree.raxmltree' outb_name];
    tree = phytreeread(filePhylo);
    tree = reroot(tree);
    

    [AMtree,WMtree,patients] = phytree2graph(tree,'_',1);
    [AMtree, WMtree,patients] = reduceTree(AMtree, WMtree,patients);
    
    patientList = sort(unique(patients));
    patientList = patientList(2:end);
    nPat = length(patientList);
    
    AMTNtrue = zeros(nPat,nPat);
    src = find(patientList == sources(i));
    AMTNtrue(src,:) = 1;
    AMTNtrue(src,src) = 0;
    sourceTrue = src;

    outb_fold = [fold filesep outb_name]; % folder with sequence data. Contains nPat files with sequences from each patient
    files = dir(fullfile(outb_fold,'*.fas'));
    divers = zeros(1,nPat);
    for j =1:size(files,1) 
       seq = fastaread([outb_fold filesep files(j).name]);
       seq = char(seq.Sequence);
       L = size(seq,2);
       tokens = strsplit(files(j).name, '_');
       patID = find(patientList == str2num(tokens{1}(3:end)));
       for col = 1:L
            [uniqueChars, ~, idx] = unique(seq(:, col));
            charCounts = accumarray(idx, 1);
            probabilities = charCounts / sum(charCounts);
            entropy = -sum(probabilities .* log2(probabilities));
            divers(patID) = divers(patID) + entropy;
       end
       divers(patID) = divers(patID)/L;
    end
    
    
    [traits,traitFreq] = pat2traits1(patients,patientList); 
    tree = digraph(AMtree);
    treeUndir = graph(AMtree+AMtree' > 0);

    predec = findPredecessors(tree,nPat);
    AMTNinferSamp = cell(1,nSamp);
    objLPSamp = -ones(1,nSamp);
    objParsSamp = -ones(1,nSamp);
    suffSamp = -ones(1,nSamp);
    parfor s = 1:nSamp
        s
        el = preferential_attachment(nPat,1);
        el = el(1:2:size(el,1),1:2);
        [ishom,map,invmap,AMTNinferCurr,objLP,objPars,suff,AMinferCurrDir,origin] = checkHomomorUnconstrMult(treeUndir,tree,cand,traits,objDeg,divers);
        AMTNinferSamp{s} = AMTNinferCurr;
        objLPSamp(s) = objLP;
        objParsSamp(s) = objPars;
        suffSamp(s) = suff;
    end
    save(outFile,'AMTNinferSamp','objLPSamp','objParsSamp','suffSamp');
end
%% 
percPars = 99;
    
for i = 1:nOut
    outb_name = outbreaks{i};
    outFile = ['results' filesep 'realHCV' filesep 'res_' outb_name '.mat'];
    load(outFile);

    filePhylo = [inputFolder filesep 'RAxML_bestTree.raxmltree' outb_name];
    tree = phytreeread(filePhylo);
    
    [AMtree,patients] = phytree2graph(tree,tree,delimeter,tokenPos);
    [AMtree, patients,~] = reduceTree(AMtree,patients);
    
    patientList = sort(unique(patients));
    patientList = patientList(2:end);
    nPat = length(patientList);
    
    AMTNtrue = zeros(nPat,nPat);
    src = find(patientList == sources(i));
    AMTNtrue(src,:) = 1;
    AMTNtrue(src,src) = 0;
    sourceTrue = src;

    p = prctile(-objSamp, percPars)
    ind = (objSamp <= -p);
    AMTNinferSamp = AMTNinferSamp(ind);
    
    AMCons = sum(cat(3, AMTNinferSamp{:}), 3)/length(AMTNinferSamp);
    AMCons(AMCons < 0.0001) = 0;
    G = graph(AMCons);
    LWidths = 5*G.Edges.Weight/max(G.Edges.Weight);
    fig = figure
    p = plot(G,'Layout','force','LineWidth',LWidths);
    p.MarkerSize = 10; 
    p.NodeFontSize = 14;
    p.NodeFontWeight = 'bold';
    ind = find((G.Edges.EndNodes(:,1) == sourceTrue) | (G.Edges.EndNodes(:,2) == sourceTrue));
    highlight(p,'Edges',ind,'EdgeColor','r');
    set(gca, 'Box', 'off');
    axis off
    xl = xlim;
    yl = ylim;
    if i == 3
        text(0.6*xl(1),0.99*yl(2),['(c)'],'FontSize',30,'FontWeight','bold','HorizontalAlignment','right','VerticalAlignment','top')
    end
    exportgraphics(fig,['figures' filesep 'consNet_' outb_name '.png'],'Resolution',600)



    [traits,traitFreq] = pat2traits1(patients,patientList); 
    [AMtree,traits,AMunLab] = contractMultLabelRelax(AMtree,traits,nPat);
    G = graph(AMunLab);
    fig = figure
    p = plot(G,'Layout','force','LineWidth',2);
    p.MarkerSize = 10; 
    p.NodeFontSize = 14;
    p.NodeFontWeight = 'bold';
    ind = find((G.Edges.EndNodes(:,1) == sourceTrue) | (G.Edges.EndNodes(:,2) == sourceTrue));
    highlight(p,'Edges',ind,'EdgeColor','r');
    set(gca, 'Box', 'off');
    axis off
    xl = xlim;
    yl = ylim;
    if i == 3
        text(0.6*xl(1),0.99*yl(2),['(b)'],'FontSize',30,'FontWeight','bold','HorizontalAlignment','right','VerticalAlignment','top')
    end
    exportgraphics(fig,['figures' filesep 'parsNet_' outb_name '.png'],'Resolution',600)
end

%% 

function [AMtreeCont,traitsCont,AMunLab] = contractMultLabelRelax(AMtree,traits,nPat)

    AMtreeUndir = (AMtree + AMtree')>0;
    AMtreeCont = [];
    traitsCont = [];

    root = find(sum(AMtree,1) == 0);
    leafs = (sum(AMtree,2)==0);
    tree = digraph(AMtree);    
    dfsorder = flip(dfsearch(tree,root))';

    labels = cell(1,numnodes(tree));
    for v = dfsorder
        if leafs(v)
            labels{v} = traits(v);
        else
            child = successors(tree,v);
            inter = intersect(labels{child(1)},labels{child(2)});
            if ~isempty(inter)
                labels{v} = inter;
            else
                labels{v} = unique(sort([labels{child(1)} labels{child(2)}]));
            end
        end
    end
    (sum(cellfun(@length, labels) == 1) - sum(leafs))/(length(AMtree)-sum(leafs))

    dfsorder = flip(dfsorder);
    for v = dfsorder
        if length(labels{v})==1
            traits(v) = labels{v};
            continue;
        end
        if (v~=root) && ismember(traits(predecessors(tree,v)),labels{v})
            traits(v) = traits(predecessors(tree,v));
        end
    end


    contract = zeros(1,length(traits));
    c = 0;
    for v = dfsorder
        if traits(v) == 0
            c = c+1;
            contract(v) = c;
        end
        if (traits(v) > 0) && (contract(v) == 0)
            AMsubtree = zeros(size(AMtree,1),size(AMtree,2));
            ind = (traits == traits(v));
            AMsubtree(ind,ind) = 1;
            AMsubtree = AMsubtree.*AMtreeUndir;
            [comps,~] = conncomp(graph(AMsubtree));
            c = c+1;
            contract(comps==comps(v)) = c;
        end
    end

    AMtreeCont = zeros(c,c);
    traitsCont = zeros(1,c);
    for i = 1:c
        for j = 1:c
            if i == j
                continue;
            end
            if sum(sum(AMtree(contract==i,contract==j))) > 0
                AMtreeCont(i,j) = 1;
            end
        end
        traitsCont(i) = traits(find(contract==i,1));
    end

    toRemove = zeros(1,length(AMtreeCont));
    for i = 1:length(AMtreeCont)
        neigh = find(AMtreeCont(i,:));
        if ~isempty(neigh)
            U = unique(traitsCont(neigh));
            for u = U
                if u == 0
                    continue;
                end
                ind = find(traitsCont(neigh) == u);
                if length(ind) > 1
                    toRemove(neigh(ind(2:end))) = 1;
                    AMtreeCont(neigh(ind(1)),:) = AMtreeCont(neigh(ind(1)),:) + sum(AMtreeCont(neigh(ind(2:end)),:),1);
                end
            end
        end
    end
    AMtreeCont = AMtreeCont(~toRemove,~toRemove);
    traitsCont = traitsCont(~toRemove);

    
    treeCont = digraph(AMtreeCont);
    root = find(sum(AMtreeCont,1) == 0);
    dfsorder = dfsearch(treeCont,root)';
    AMunLab = zeros(nPat,nPat);
    for v = dfsorder
        if v ~= root
            p = predecessors(treeCont,v);
            if (traitsCont(p) > 0) && (traitsCont(v) > 0)
                AMunLab(traitsCont(p),traitsCont(v)) = 1;
                AMunLab(traitsCont(v),traitsCont(p)) = 1;
            end
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

%% 
