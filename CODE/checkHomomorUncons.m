function [ishom,AM,obj,origin] = checkHomomorUncons(tree,cand,traits,divers,timeLimit)

eps = 0.001;
n = numnodes(cand);

degCand = degree(cand);
degTree = zeros(1,length(traits));
if isempty(divers)
    for i = 1:length(traits)
        node_i = find(traits==i);
        for j = 1:length(node_i)
            neigh = successors(tree,node_i(j));
            degTree(i) = degTree(i) + length(neigh);
        end
    end
    degTree = degTree(traits > 0);
else
    degTree = divers;
end

Atree = adjacency(tree);
treeUndir = graph(Atree + Atree' > 0);
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
obj = Inf;


id = randi(1000000);
filename = ['ILP_' int2str(id) '.lp'];
fileID = fopen(filename,'w');


    fprintf(fileID,'%s\n','maximize');    
    for i = 1:n
        for j = 1:n
            Xij = ['X' int2str(i) ',' int2str(j)];
            if sum(degTree) ~= 0
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
            obj = Inf;
            AM = [];
            return;
        end
    end
end
labels = zeros(1,numnodes(tree));
dfsorder = flip(dfsorder);
obj = 0;
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
        obj = obj + 1;
    end
end
