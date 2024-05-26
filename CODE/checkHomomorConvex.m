function [ishom,AMinfer,obj,origin] = checkHomomorConvex(tree,cand,traits,toOptimizeComp,divers)

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
    quasidegTree(i) = outdegree(tree,i);
end
quasidegTree = quasidegTree(indTraits);

n = length(degCand);
N = numnodes(tree);
r_tree = find(indegree(tree)==0);
predec_tree = findPredecessors(tree,n);


if sum(degCand==1) == n-1
    ishom = true;
    if ~isempty(divers)
        centr = find(divers==max(divers),1);
    else
        centr = randi(n);
    end
    AMinfer = zeros(n,n);
    AMinfer(centr,:) = 1;
    AMinfer(:,centr) = 1;
    AMinfer(centr,centr) = 0;
    obj = quasidegTree(centr)*(n-1);
    origin = centr;
    return;
end

neighborhoods = arrayfun(@(x) find(AM_cand(x,:)), 1:n, 'UniformOutput', false);
neighborhoods_str = cellfun(@mat2str, neighborhoods, 'UniformOutput', false);
[~,ia_cand,ic_cand] = unique(neighborhoods_str,'stable');
n = length(ia_cand);
counts_cand = histcounts(ic_cand, 1:n+1);
AM_cand_cont = AM_cand(ia_cand,ia_cand);
cand_cont = graph(AM_cand_cont);
degCandCont = degree(cand_cont);
m = numedges(cand_cont);

AM_tree_aux = [AM_tree AM_tree'];
neighborhoods = arrayfun(@(x) find(AM_tree_aux(x,:)), 1:N, 'UniformOutput', false);
neighborhoods_str = cellfun(@mat2str, neighborhoods, 'UniformOutput', false);
[~,ia_tree,ic_tree] = unique(neighborhoods_str,'stable');
N = length(ia_tree);
counts_tree = histcounts(ic_tree, 1:N+1);
AM_tree_cont = AM_tree(ia_tree,ia_tree);
tree_cont = digraph(AM_tree_cont);
traits_cont = traits(ia_tree);

r_tree_cont = find(indegree(tree_cont)==0);
orderDown_cont = (dfsearch(tree_cont,r_tree_cont))';
orderUp = flip(orderDown_cont);
outdeg_tree_cont = outdegree(tree_cont);
AM_cand_cont_ref = AM_cand_cont + eye(n,n);



    possHomomor = cell(1,N); %partial homomorphisms at each node of phylogeny in the form (w,X,count), where X represented by a characteristic vector
    generateHomomor = cell(1,N); % ids of child homomorphisms forming each parent homomorphism
    coveredHomomor = cell(1,N); % vertices assigned to labelled nodes in partial homomorphisms

    nPartialDesc = zeros(n,n);
    E = cand_cont.Edges.EndNodes;
    for e = 1:m
        A = AM_cand_cont;
        A(E(e,1),E(e,2)) = 0;
        A(E(e,2),E(e,1)) = 0;
        H = graph(A);
        [comps,~] = conncomp(H);
        nPartialDesc(E(e,1),E(e,2)) = sum(counts_cand( comps==comps(E(e,2) )));   % #descendants of E(e,2) | E(e,1) is a parent
        nPartialDesc(E(e,2),E(e,1)) = sum(counts_cand( comps==comps(E(e,1) ))); 
    end
    lf_attach = zeros(1,n);
    for i = 1:n
        if degree(cand_cont,i) == 1
            lf_attach(neighbors(cand_cont,i)) = i;
        end
    end
    intern_cand_cont = (find(degCandCont > 1))';

    for a = orderUp
        if outdeg_tree_cont(a) == 0
            if counts_tree(a) == 1
                homs = zeros(n,n+1);
                homs(:,1) = (1:n)';
                homs(:,2:end) = eye(n,n);
            else
                p = predecessors(tree_cont,a);
                possImagesLeafs  = find((counts_cand>=counts_tree(a)) & (degCandCont' == 1));
                if traits_cont(p) == 0
                    possImagesIntern = find(lf_attach>0);
                    possImagesInternLeafs = lf_attach(possImagesIntern);
                    possImagesIntern = possImagesIntern(counts_cand(possImagesInternLeafs)>=counts_tree(a)-1);
                else
                    possImagesIntern = [];
                end
    
                homs = zeros(length(possImagesLeafs)+length(possImagesIntern),n+1);
                homs(:,1) = [possImagesLeafs'; possImagesIntern'];
                homs(sub2ind(size(homs), 1:length(possImagesLeafs), possImagesLeafs+1)) = counts_tree(a);
                if ~isempty(possImagesIntern)
                    homs(sub2ind(size(homs), (length(possImagesLeafs)+1):(length(possImagesLeafs)+length(possImagesIntern)),... 
                        possImagesIntern+1)) = 1;
                    homs(sub2ind(size(homs), (length(possImagesLeafs)+1):(length(possImagesLeafs)+length(possImagesIntern)),... 
                        possImagesInternLeafs+1)) = counts_tree(a)-1;
                end
            end


            cov = homs(:,2:end);            
            possHomomor{a} = homs;
            coveredHomomor{a} = cov;

            continue;
        end

        child = successors(tree_cont,a);
        leaf_tree = child(outdeg_tree_cont(child)==0);
        leaf_tree_id = find(outdeg_tree_cont(child)==0);
        possHomomor_child = possHomomor(child);
        coveredHomomor_child = coveredHomomor(child);
        poss_hom_sizes = cellfun(@(x) size(x, 1), possHomomor_child);
        homs_a = zeros(prod(poss_hom_sizes),n+1);
        generate_homs_a = zeros(prod(poss_hom_sizes),length(child));
        assignLeaf_a = zeros(prod(poss_hom_sizes),n);
        viol_a = zeros(prod(poss_hom_sizes),1);
        
        homRootsChild = cellfun(@(x) x(:,1),possHomomor_child,'UniformOutput', false);
        blockPos = [0 cumsum(poss_hom_sizes)];
        homConcat= vertcat(possHomomor_child{:});
        coveredHomomorConcat = vertcat(coveredHomomor_child{:});
        childIDConcat = repelem(1:length(child), poss_hom_sizes);

        countNewH = 0;
        if a~=r_tree_cont
            poss_w = intern_cand_cont;
        else
            poss_w = 1:n;
        end
        for w = poss_w

            allnonadj_w = cellfun(@(x) all(AM_cand_cont_ref(w,x) == 0), homRootsChild);
            if any(allnonadj_w)
                continue;
            end

            if lf_attach(w)>0
                cap_leaf_w = counts_cand(lf_attach(w));
                lf_w = lf_attach(w);
            else
                cap_leaf_w = 0;
                lf_w = [];
            end
 
            if ~isempty(leaf_tree)
                ind_adj_noleaf = (childIDConcat~=leaf_tree_id) & AM_cand_cont_ref(w,homConcat(:,1));
                ind_adj_leaf = (childIDConcat==leaf_tree_id) & ... 
                    ((coveredHomomorConcat*AM_cand_cont_ref(:,w)==counts_tree(leaf_tree)))';
            else
                 ind_adj_noleaf = AM_cand_cont_ref(w,homConcat(:,1))>0;
                 ind_adj_leaf = false(1,length(childIDConcat));
            end
            ind_subtree_cov = (homConcat(:,1)~=w) & ((homConcat(:,1)==lf_attach(w)) |  ...
                (sum(coveredHomomorConcat,2)==(nPartialDesc(w,homConcat(:,1)))'));
            if traits_cont(a) == 0
                hom_compat_w = find((ind_adj_leaf | ind_adj_noleaf) & (ind_subtree_cov' | (homConcat(:,1)==w)'));
            else
                hom_compat_w = find((ind_adj_leaf | ind_adj_noleaf) & (ind_subtree_cov' | (homConcat(:,1)==w)')  & ... 
                    (~coveredHomomorConcat(:,w))');
            end

            homConcat_w = homConcat(hom_compat_w,:);
            childIDConcat_w = childIDConcat(hom_compat_w);

            cnt = zeros(1,outdeg_tree_cont(a));
            cnt(childIDConcat_w) = 1;
            if sum(cnt) < outdeg_tree_cont(a)
                continue;
            end

            homConcat_w(:,w+1) = 0;
            homConcat_w(:,lf_w+1) = 0;
            M_intersect_w = (homConcat_w(:,2:end)*(homConcat_w(:,2:end))');
            for i = 1:length(child)
                ind = childIDConcat_w==i;
                M_intersect_w(ind,ind) = 1;
            end
            homConcat_w = homConcat(hom_compat_w,:);
            if ~isempty(lf_w)
                hom_rootwl = find(homConcat_w(:,1)==w | homConcat_w(:,1)==lf_w);
            else
                hom_rootwl = find(homConcat_w(:,1)==w);
            end
            if length(hom_rootwl) > 1
                pairs = nchoosek(hom_rootwl,2);
                if ~isempty(lf_w)
                    ind = (homConcat_w(pairs(:,1),lf_w+1)+homConcat_w(pairs(:,2),lf_w+1)  > counts_cand(lf_w));
                    pairs = pairs(ind,:);
                end
                M_intersect_w(sub2ind(size(M_intersect_w), pairs(:,1), pairs(:,2))) = 1;
                M_intersect_w(sub2ind(size(M_intersect_w), pairs(:,2), pairs(:,1))) = 1;
            end             
            AM_w = (M_intersect_w == 0);

            trans = getTransversalsClique(AM_w,length(child));
            for i = 1:length(trans)%size(trans,2)
                tr = hom_compat_w(trans{i});
                hom_join = homConcat(tr,:);

                copies_w = tr(homConcat(tr,1) == w);
                cnt_leaf_w = max([0 sum(homConcat(copies_w,lf_w+1))]);

                if cnt_leaf_w > cap_leaf_w
                    continue;
                end

                cov = sum(coveredHomomorConcat(tr,:),1);
                if traits_cont(a) > 0
                    cov(w) = 1;
                end
                hm = sum(hom_join(:,2:end),1);
                hm(w) = 1;

                ind = true(1,n);
                ind([w lf_w]) = 0;
                if any(hm(ind) ~= cov(ind))
                    continue;
                end

                ind = (cov>0);
                ind(lf_w) = 0; 
                if any(counts_cand(ind)~= cov(ind)) 
                    continue;
                end
                countNewH = countNewH + 1;
                homs_a(countNewH,2:end) = hm;
                homs_a(countNewH,1) = w;
                generate_homs_a(countNewH,:) = tr - blockPos(childIDConcat(tr));
                assignLeaf_a(countNewH,:) = cov;
                viol_a(countNewH) = sum(hom_join(:,1)==w);
            end
        end

        if countNewH == 0
            return;
        end
        homs_a = homs_a(1:countNewH, :);
        generate_homs_a = generate_homs_a(1:countNewH,:);
        assignLeaf_a = assignLeaf_a(1:countNewH,:);
        viol_a = viol_a(1:countNewH);

        [~, ~, ic] = unique(homs_a, 'rows', 'stable');
        minViolValues = accumarray(ic, viol_a, [], @min);
        ind = arrayfun(@(i) find((viol_a == minViolValues(i))&(ic==i), 1, 'first'), 1:length(minViolValues));        
        homs_a = homs_a(ind, :);
        generate_homs_a = generate_homs_a(ind,:);
        assignLeaf_a = assignLeaf_a(ind,:);
     
        possHomomor{a} = homs_a;
        generateHomomor{a} = generate_homs_a;
        coveredHomomor{a} = assignLeaf_a;

    end

    if isempty(possHomomor{r_tree_cont})
        return;
    end


    nSol = size(possHomomor{r_tree_cont},1);
    for s = 1:nSol
        leafRoot = 0;
        homomor_select = zeros(1,N);
        mapping_cont = zeros(1,N);
        for a = orderDown_cont
            child = successors(tree_cont,a);            
            homs_a = possHomomor{a};
            if a == r_tree_cont
                homomor_select(a) = s;
            end
            hom = homs_a(homomor_select(a),:);
            mapping_cont(a) = hom(1);
            if isempty(child)
                continue;
            end
            if (a == r_tree_cont) && (degCandCont(hom(1))==1) && (length(child) > 1)
                leafRoot = hom(1);
            end
            homomor_select(child) = generateHomomor{a}(homomor_select(a),:);
        end

        
        N_full = numnodes(tree);
        mapping_cont1 = zeros(1,N_full);
        for j = 1:length(mapping_cont)
            mapping_cont1(ic_tree == j) = mapping_cont(j);
        end

        nodes2 = (find((outdegree(tree)==2) & (traits == 0)))';
        for v = nodes2
            child = successors(tree,v);
            if sum(outdegree(tree,child)) > 0
                continue;
            end
            if (mapping_cont1(v) == mapping_cont1(child(1))) && (mapping_cont1(v) == mapping_cont1(child(2)))
                hm = possHomomor{ic_tree(v)}(homomor_select(ic_tree(v)),:);
                hm(mapping_cont1(v)+1) = 0;
                img = find(hm(2:end));
                if ~isempty(divers)
                    if divers(traits(child(1))) > divers(traits(child(2)))
                        mapping_cont1(child(2)) = img;
                    else
                        mapping_cont1(child(1)) = img;
                    end
                else
                    ii = randi(2);
                    mapping_cont1(child(ii)) = img;
                end
            end
        end
    
        mapping = zeros(1,N_full);
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
    
        n_full = numnodes(cand);
        AMinfer_curr = zeros(n_full,n_full);
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
        if toOptimizeComp
            obj_new = obj_new/(n^2) - checkSamplViolate(tree,AMinfer_curr,traits,predec_tree);
        end
        if obj_new > obj
            obj = obj_new;
            AMinfer = AMinfer_curr;
            origin = origin_curr;
        end
    end


