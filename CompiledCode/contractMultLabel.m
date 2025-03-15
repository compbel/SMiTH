function [AMtreeCont,traitsCont] = contractMultLabel(AMtree,traits,nPat)

    root = find(sum(AMtree,1) == 0);
    tree = digraph(AMtree);    
    dfsorder = flip(dfsearch(tree,root))';

    recCommAnc = zeros(1,nPat);
    nDescTrait = zeros(length(dfsorder),nPat);
    for v = dfsorder
        if traits(v) ~= 0
            nDescTrait(v,traits(v)) = 1;
        else
            child = successors(tree,v);
            nDescTrait(v,:) = sum(nDescTrait(child,:),1);
        end
    end    
    for tr = 1:nPat
        nTrLeaves = sum(traits == tr,2);
        ancTr = nDescTrait(dfsorder,tr);
        recCommAnc(tr) = dfsorder(find(ancTr == nTrLeaves,1,"first"));
    end
    for i = 1:nPat
        trLeaves = find(traits == i);
        for v = trLeaves
            P = shortestpath(tree,recCommAnc(i),v);
            for w = P
                if (traits(w) > 0) && (traits(w) ~= i)
                    ['No solution']
                    return;
                else
                    traits(w) = i;
                end
            end
        end
    end

    vertCon = true(1,length(traits));
    for tr = 1:nPat
        vertLabel = (traits == tr);
        if sum(vertLabel) >= 2
            outneigh = find((traits~=tr)&(sum(AMtree(vertLabel,:),1)>0));
            AMtree(recCommAnc(tr),outneigh) = 1;
        end
        vertLabel(recCommAnc(tr)) = 0;
        vertCon(vertLabel) = 0;
    end
    AMtreeCont = AMtree(vertCon,vertCon);
    traitsCont = traits(vertCon);
