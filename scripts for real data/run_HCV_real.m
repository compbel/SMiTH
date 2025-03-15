clear;

inputFolder = 'HCV_real';
outbreaks = {'AI','AW'};
sources = [4 2];
nOut = length(outbreaks);
nSamp = 100000;
timeLimit = 600;

delimeter = '_';
tokenPos = 1;
sampGenerator = @randTreePrefAttach;
constr = 'unconstrained'; 

%% 
fscoresCons = zeros(1,length(outbreaks));
fscoresTotal = zeros(1,length(outbreaks));
for i = 1:nOut
    outb_name = outbreaks{i};
    outFile = ['results' filesep 'realHCV' filesep 'res_' outb_name '.mat'];
    filePhylo = [inputFolder filesep 'RAxML_bestTree.raxmltree' outb_name];
    fileSeq = [inputFolder filesep  outbreaks{i} '_all.fas'];


    [migrSamp,objSamp,originSamp,consensus, siteList] = migrationSampler(filePhylo,sampGenerator,...
    nSamp,constr,timeLimit,fileSeq,delimeter,tokenPos);
    save(outFile,'migrSamp','objSamp','suffSamp');

    [];
end
%% 
% calculate results
percLP = 99;
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

    AMTNinferSamp = cellfun(@adjacency, migrSamp, 'UniformOutput', false);

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