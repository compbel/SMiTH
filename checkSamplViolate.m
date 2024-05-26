function numViol = checkSamplViolate(tree,AMcand,traits,predec)

numViol = 0;
intern = (find(outdegree(tree)>0))';
for p = intern
    if traits(p) > 0
        patsConstrBegin = p;
    else
        patsConstrBegin = (find(predec(:,p)>0))';
    end
    for i = patsConstrBegin
        child = (successors(tree,p))';
        for c = child
            if predec(i,c) ~= 1
                patsConstrEnd = (find(predec(:,c)))';
                if length(patsConstrEnd) <= 1
                    continue;
                end
                numAdjSubtree = 0;
                for j = patsConstrEnd
                    if AMcand(i,j) == 1
                        numAdjSubtree = numAdjSubtree + 1;
                    end
                end
                if numAdjSubtree >= 2
                    numViol = numViol + 1;
                end
            end
        end
    end
end


