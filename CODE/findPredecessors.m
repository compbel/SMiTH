function predec = findPredecessors(tree,nPat)

nNodes = numnodes(tree);
predec = (adjacency(transclosure(tree)))'+ eye(nNodes,nNodes);
predec = predec(1:nPat,:);
