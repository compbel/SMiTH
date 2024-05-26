function [ishom,AMinfer,obj,origin] = checkHomomorCompact(tree,cand,traits)

ishom = false;
AMinfer = [];
obj = -Inf;
origin = [];

indTraits = find(traits > 0);
degCand = degree(cand);
AM_cand = adjacency(cand);
AM_tree = adjacency(tree);
quasidegTree = zeros(1,length(traits));
for i = 1:length(traits)
    neigh = successors(tree,i);
    % degTree(i) = sum(indTraits(neigh));
    quasidegTree(i) = length(neigh);
end
quasidegTree = quasidegTree(indTraits);

n = length(degCand);
N = numnodes(tree);
r_tree = find(indegree(tree)==0);

if sum(degCand==1) == n-1
    intern = outdegree(tree,1:N)>0;
    intern_neigh = zeros(1,N);
    for i = 1:N
        intern_neigh(i) = sum(intern(successors(tree,i)));
    end
    if max(intern_neigh(intern)) <= 1
        u = find(intern + intern_neigh == 1);
        if length(u) == 1
            ishom = true;
            AMinfer = zeros(n,n);
            AMinfer(u,:) = 1;
            AMinfer(:,u) = 1;
            AMinfer(u,u) = 0;
            obj = quasidegTree(u)*(n-1);
            origin = u;
            return;
        end
    end
end

neighborhoods = arrayfun(@(x) find(AM_cand(x,:)), 1:n, 'UniformOutput', false);
neighborhoods_str = cellfun(@mat2str, neighborhoods, 'UniformOutput', false);
[~,ia_cand,ic_cand] = unique(neighborhoods_str,'stable');
n = length(ia_cand);
counts_cand = histcounts(ic_cand, 1:n+1);
% AM_cand_cont = AM_cand(ia,ia).*(counts'*counts);
AM_cand_cont = AM_cand(ia_cand,ia_cand);
cand_cont = graph(AM_cand_cont);
degCandCont = degree(cand_cont);
m = numedges(cand_cont);
neighborhoods_leaf = arrayfun(@(x) find(AM_cand_cont(x,:)&(degCandCont'==1)), 1:n, 'UniformOutput', false);

AM_tree_undir = (AM_tree + AM_tree') > 0;
neighborhoods = arrayfun(@(x) find(AM_tree_undir(x,:)), 1:N, 'UniformOutput', false);
neighborhoods_str = cellfun(@mat2str, neighborhoods, 'UniformOutput', false);
for i = 1:N
    p = find(AM_tree_undir(:,i));
    if (length(p) == 1) && (traits(p) == 0)
        neighborhoods_str{i} = mat2str([p,i]);
    end
end
[~,ia_tree,ic_tree] = unique(neighborhoods_str,'stable');
N = length(ia_tree);
counts_tree = histcounts(ic_tree, 1:N+1);
AM_tree_cont = AM_tree(ia_tree,ia_tree);
tree_cont = digraph(AM_tree_cont);
traits_cont = traits(ia_tree);
predec_tree_cont = findPredecessors(tree_cont,N);
predec_tree_cont = diag(counts_tree)*predec_tree_cont;
predec_tree_cont = predec_tree_cont(traits_cont>0,:);

r_tree_cont = find(indegree(tree_cont)==0);
orderDown_cont = (dfsearch(tree_cont,r_tree_cont))';
orderUp_cont = flip(orderDown_cont);
outdeg_tree_cont = outdegree(tree_cont);



    possHomomor = cell(1,N); %partial homomorphisms at each node of phylogeny
    transHomomor = cell(1,N); % transversals of childern of each node that create the corresponding partial homomorphism

    nPartialDesc = zeros(n,n);
    E = cand_cont.Edges.EndNodes;
    for e = 1:m
        H = rmedge(cand_cont,E(e,1),E(e,2));
        [comps,compsizes] = conncomp(H);
        nPartialDesc(E(e,1),E(e,2)) = sum(counts_cand( comps==comps(E(e,2) )));   % #descendants of E(e,2) | E(e,1) is a parent
        nPartialDesc(E(e,2),E(e,1)) = sum(counts_cand( comps==comps(E(e,1) ))); 
    end
    lf_attach = zeros(1,n);
    for i = 1:n
        if degree(cand_cont,i) == 1
            lf_attach(neighbors(cand_cont,i)) = i;
        end
    end

    for a = orderUp_cont
        % if a == 3
        %     [];
        % end
        if outdeg_tree_cont(a) == 0
            map = containers.Map('KeyType', 'char', 'ValueType', 'any');
            if counts_tree(a) > 1
                possImages  = find(counts_cand>=counts_tree(a));
                cl = counts_tree(a);
            else
                possImages = 1:n;
                cl = 0;
            end
            cl = 0;
            for i = 1:length(possImages)
                map(jsonencode([possImages(i) cl])) = [possImages(i) cl];
            end
            possHomomor{a} = map;
            continue;
        end
        child = successors(tree_cont,a);
        homs_a = containers.Map('KeyType', 'char', 'ValueType', 'any');
        transversals_a = containers.Map('KeyType', 'char', 'ValueType', 'any');
        if traits_cont(a) == 0
            hom1 = possHomomor{child(1)};
            hom2 = possHomomor{child(2)};

            keys1 = keys(hom1);
            keys2 = keys(hom2);
            for i = 1:length(keys1)
                for j = 1:length(keys2)
                    % h: root|child1,...,child_n|leaf
                    h1 = hom1(keys1{i});
                    h2 = hom2(keys2{j});
                    cl1 = h1(end);
                    cl2 = h2(end);
                    h1 = h1(1:end-1);
                    h2 = h2(1:end-1);
                    hnew = {};
                    transnew = {};
                    if (AM_cand_cont(h1(1),h2(1)) == 1) && (degCandCont(h2(1)) == 1) && (cl1 < counts_cand(h2(1))) && (~ismember(h1(1),h2(2:end)))
                        if ismember(h2(1),h1(2:end))
                            hnew = [hnew [h1 cl1+1]];
                        else
                            hnew = [hnew [h1(1) sort([h1(2:end) h2(1)]) cl1+1] ];
                        end
                        transnew = [transnew [h1(1), h2(1)]];
                        if (a == r_tree_cont) && (cl1 == counts_cand(h2(1))-1)
                            hnew = [hnew [h2(1) h1(1) 0]];
                            transnew = [transnew [h1(1), h2(1)]];
                        end
                    end
                    if (AM_cand_cont(h2(1),h1(1)) == 1) && (degCandCont(h1(1)) == 1) && (cl2 < counts_cand(h1(1))) && (~ismember(h2(1),h1(2:end)))
                        if ismember(h1(1),h2(2:end))
                            hnew = [hnew [h2 cl2+1]];
                        else
                            hnew = [hnew [h2(1) sort([h2(2:end) h1(1)]) cl2+1] ];
                        end
                        transnew = [transnew [h1(1),h2(1)] ];
                        if (a == r_tree_cont) && (cl2 == counts_cand(h1(1))-1)
                            hnew = [hnew [h1(1) h2(1) 0]];
                            transnew = [transnew [h1(1), h2(1)]];
                        end
                    end
                    if (AM_cand_cont(h1(1),h2(1)) == 1) && (degCandCont(h1(1)) > 1) && (degCandCont(h2(1)) > 1) && (~ismember(h2(1),h1(2:end))) && (~ismember(h1(1),h2(2:end)))
                        h1_noleaf = setdiff(h1(2:end),neighborhoods_leaf{h1(1)},'stable');
                        if sum(predec_tree_cont(:,a)) == sum(nPartialDesc(h1(1),h1_noleaf)) + cl1 + nPartialDesc(h1(1),h2(1)) + 1
                            hnew = [hnew [h1(1) sort([h1(2:end) h2(1)]) cl1]];
                            transnew = [transnew [h1(1), h2(1)] ];
                        else
                            [];
                        end
                    end
                    if (AM_cand_cont(h2(1),h1(1)) == 1) && (degCandCont(h1(1)) > 1) && (degCandCont(h2(1)) > 1) && (~ismember(h1(1),h2(2:end))) && (~ismember(h2(1),h1(2:end)))
                        h2_noleaf = setdiff(h2(2:end),neighborhoods_leaf{h2(1)},'stable');
                        if sum(predec_tree_cont(:,a)) == sum(nPartialDesc(h2(1),h2_noleaf)) + cl2 + nPartialDesc(h2(1),h1(1)) + 1
                            hnew = [hnew [h2(1) sort([h2(2:end) h1(1)]) cl2]];
                            transnew = [transnew [h1(1),h2(1)] ];
                        else
                            [];
                        end
                    end
                    if ~isempty(hnew)
                        for k = 1:length(hnew)
                            keynew = jsonencode(hnew{k});
                            if ~isKey(homs_a,keynew)
                                homs_a(keynew) = hnew{k};
                                sources = containers.Map('KeyType', 'char', 'ValueType', 'any');
                                sources(jsonencode(transnew{k})) = transnew{k};
                                transversals_a(keynew) = sources;
                            else
                                trans_prev = transversals_a(keynew);
                                trans_prev(jsonencode(transnew{k})) = transnew{k};
                                transversals_a(keynew) = trans_prev;
                            end
                        end
                    end
                end
            end
            if isempty(keys(homs_a))
                break;
            end
            possHomomor{a} = homs_a;
            transHomomor{a} = transversals_a;
        else
            if (a == r_tree_cont) && (length(child) == 1)
                hom = possHomomor{child};
                keysc = keys(hom);
                for j = 1:length(keysc)
                    h = hom(keysc{j});
                    nb = neighbors(cand_cont,h(1));
                    lf = nb(degCandCont(nb)==1);
                    if (~isempty(lf)) && (h(end)==counts_cand(lf)-1) && (sum(predec_tree_cont(:,a))==sum(nPartialDesc(h(1),:))+1)
                        hnew = [lf h(1) 0];
                        keynew = jsonencode(hnew);
                        homs_a(keynew) = hnew;
                        sources = containers.Map('KeyType', 'char', 'ValueType', 'any');
                        sources(jsonencode(h(1))) = h(1);
                        transversals_a(keynew) = sources;
                    end
                end
                if isempty(keys(homs_a))
                    return;
                end
                possHomomor{a} = homs_a;
                transHomomor{a} = transversals_a;
                continue;
            end
            sets = cell(1,length(child));
            leaf_tree = child(outdeg_tree_cont(child)==0);
            for i = 1:length(child)
                b = child(i);
                sb = zeros(1,n);
                homb = possHomomor{b};
                keysb = keys(possHomomor{b});
                for j = 1:length(keysb)
                    hb = homb(keysb{j});
                    if (length(child)>1) && (a ~= r_tree_cont) && (lf_attach(hb(1)) > 0) && (hb(end) < counts_cand(lf_attach(hb(1))))
                        continue;
                    end
                    sb(hb(1)) = 1;
                end
                sets{i} = find(sb);
            end
            trans = getTransversals(sets,n);
            for i = 1:size(trans,2)
                tr = (trans(:,i))';
                deg_trans = sum(AM_cand_cont(tr,:),1);
                W = find(deg_trans == length(child));
                if isempty(W)
                    continue;
                end
                for j = 1:length(W) 
                    w = W(j);
                    nonLeaf_cand_cont = tr(degCandCont(tr)>1);
                    if isempty(leaf_tree)
                        cl = 0;
                    else
                        cl = counts_tree(leaf_tree);
                    end
                    if sum(predec_tree_cont(:,a)) ~= sum(nPartialDesc(w,nonLeaf_cand_cont)) + cl + 1
                        continue;
                    end
                
                    hnew = [w sort(tr) cl];
                    keynew = jsonencode(hnew);
                    if ~isKey(homs_a,keynew)
                        homs_a(keynew) = hnew;
                        sources = containers.Map('KeyType', 'char', 'ValueType', 'any');
                        sources(jsonencode(tr)) = tr;
                        transversals_a(keynew) = sources;
                    else
                        trans_prev = transversals_a(keynew);
                        trans_prev(jsonencode(tr)) = tr;
                        transversals_a(keynew) = trans_prev;
                    end
                end
            end
            if isempty(keys(homs_a))
                return;
            end
            possHomomor{a} = homs_a;
            transHomomor{a} = transversals_a;
        end
    end

    if isempty(possHomomor{r_tree_cont})
        return;
    end

    q = 1;
    nSol = possHomomor{r_tree_cont}.Count;
    for s = 1:nSol
        leafRoot = 0;
        keys_homomor_select = cell(1,N);
        mapping_cont = zeros(1,N);
        for a = orderDown_cont
            child = successors(tree_cont,a);            
            homs_a = possHomomor{a};
            trans_a = transHomomor{a};
            if a == r_tree_cont
                keys_a = keys(homs_a);
                keys_homomor_select{a} = keys_a{s};
            end
            hom = homs_a(keys_homomor_select{a});
            cl = hom(end);
            hom = hom(1:end-1);
            mapping_cont(a) = hom(1);
            if isempty(child)
                continue;
            end
            tr = trans_a(keys_homomor_select{a});
            keys_tr = keys(tr);
            tr = tr(keys_tr{q});
            if (a == r_tree_cont) && (degCandCont(hom(1))==1)
                if length(child) > 1
                    leafRoot = hom(1);
                    if tr(1) == hom(1)
                        hom_prev1 = [setdiff(hom,tr(2),'stable') cl];
                        if counts_cand(hom(1)) > 1
                            hom_prev2 = [tr(2) (sort( neighbors(cand_cont,tr(2)) ))' counts_cand(hom(1))-1];
                        else
                            hom_prev2 = [tr(2) (sort( setdiff(neighbors(cand_cont,tr(2)),tr(1),'stable') ))' 0];
                        end
                    end
                    if tr(2) == hom(1)
                        hom_prev2 = [setdiff(hom,tr(1),'stable') cl];
                        if counts_cand(hom(1)) > 1
                            hom_prev1 = [tr(1) (sort( neighbors(cand_cont,tr(1)) ))' counts_cand(hom(1))-1];
                        else
                            hom_prev1 = [tr(1) (sort( setdiff(neighbors(cand_cont,tr(1)),tr(2),'stable') ))' 0];
                        end
                    end
                    keys_homomor_select{child(1)} = jsonencode(hom_prev1);
                    keys_homomor_select{child(2)} = jsonencode(hom_prev2);
                else
                    nb = (neighbors(cand_cont,tr(1)))';
                    if counts_cand(hom(1)) > 1
                        hom_prev = [tr(1) nb counts_cand(hom(1))-1];
                    else
                        hom_prev = [tr(1) setdiff(nb,hom(1)) 0];
                    end
                    keys_homomor_select{child} = jsonencode(hom_prev);
                end
                continue;
            end
            if traits_cont(a) == 0
                if tr(1) == hom(1)
                    if degCandCont(tr(2)) == 1
                        if cl == 1
                            hom_prev1 = [setdiff(hom,tr(2),'stable') 0];
                        else
                            hom_prev1 = [hom cl-1];
                        end
                        hom_prev2 = [tr(2) 0];
                    else
                        hom_prev1 = [setdiff(hom,tr(2),'stable') cl];
                        hom_prev2 = [tr(2) (sort( setdiff(neighbors(cand_cont,tr(2)),tr(1),'stable') ))'];
                        lf = neighborhoods_leaf{tr(2)};
                        if ~isempty(lf)
                            hom_prev2 = [hom_prev2 counts_cand(lf)];
                        else
                            hom_prev2 = [hom_prev2 0];
                        end
                    end
                end
                if tr(2) == hom(1)
                    if degCandCont(tr(1)) == 1
                        if cl == 1
                            hom_prev2 = [setdiff(hom,tr(1),'stable') 0];
                        else
                            hom_prev2 = [hom cl-1];
                        end
                        hom_prev1 = [tr(1) 0];
                    else
                        hom_prev2 = [setdiff(hom,tr(1),'stable') cl];
                        hom_prev1 = [tr(1) (sort( setdiff(neighbors(cand_cont,tr(1)),tr(2),'stable') ))'];
                        lf = neighborhoods_leaf{tr(1)};
                        if ~isempty(lf)
                            hom_prev1 = [hom_prev1 counts_cand(lf)];
                        else
                            hom_prev1 = [hom_prev1 0];
                        end
                    end
                end
                keys_homomor_select{child(1)} = jsonencode(hom_prev1);
                keys_homomor_select{child(2)} = jsonencode(hom_prev2);
            else
                for j = 1:length(child)
                    if outdeg_tree_cont(child(j)) > 0
                        hom_prev = [tr(j) sort( setdiff(neighbors(cand_cont,tr(j)),hom(1),'stable') )'];
                        lf_cand = neighborhoods_leaf{tr(j)};
                        if ~isempty(lf_cand)
                            hom_prev = [hom_prev counts_cand(lf_cand)];
                        else
                            hom_prev = [hom_prev 0];
                        end
                    else
                        hom_prev = [tr(j) 0];
                    end
                    keys_homomor_select{child(j)} = jsonencode(hom_prev);
                end
            end
        end
        
        N = numnodes(tree);
        mapping_cont1 = zeros(1,N);
        for j = 1:length(mapping_cont)
            mapping_cont1(ic_tree == j) = mapping_cont(j);
        end
    
        mapping = zeros(1,N);
        for i = 1:n
            allTwins_i = find(ic_cand == i);
            ind = find(mapping_cont1 == i);
            if i == leafRoot
                ind = setdiff(ind,r_tree);
                child = successors(tree,r_tree);
                leaf_neigh = child(outdegree(tree,child) == 0);
            end
            mapping(ind) = allTwins_i;
            if i == leafRoot
                mapping(r_tree) = mapping(leaf_neigh);
            end
        end
    
        n = numnodes(cand);
        AMinfer_curr = zeros(n,n);
        orderDown = (dfsearch(tree,r_tree))';
        for a = orderDown
            if a == r_tree
                origin_curr = mapping(a);
                continue;
            end
            p = predecessors(tree,a);
            if (mapping(a)~=mapping(p))
                i = find(mapping(indTraits)==mapping(p));
                j = find(mapping(indTraits)==mapping(a));
                AMinfer_curr(indTraits(i),indTraits(j)) = 1;
                AMinfer_curr(indTraits(j),indTraits(i)) = 1;
            end
        end
        ishom = true;
        obj_new = quasidegTree(indTraits)*degCand(mapping(indTraits));
        if obj_new > obj
            obj = obj_new;
            AMinfer = AMinfer_curr;
            origin = origin_curr;
        end
    end

