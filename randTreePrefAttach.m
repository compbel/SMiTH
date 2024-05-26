function G = randTreePrefAttach(n)
    edges = [1 2];
    deg = [1 1];    
    for i = 3:n
        prob = deg / sum(deg);
        cumProb = cumsum(prob);        
        r = rand;
        j = find(cumProb >= r, 1);        
        edges = [edges; i j];
        deg(i) = 1;
        deg(j) = deg(j) + 1;
    end    
    G = graph(edges(:,1), edges(:,2));
end