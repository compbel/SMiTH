function cliques = getTransversalsClique(AM,k)

cliques = maximalCliques(AM);
ind = cellfun(@(x) length(x)==k,cliques);
cliques = cliques(ind);
[];

