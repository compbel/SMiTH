function trans = getTransversals(sets,n)
k = length(sets);
bipart = zeros(n,k);
for i = 1:k
    bipart(sets{i},i) = 1;
end
AM = [zeros(k,k) bipart'; bipart zeros(n,n)];
G = graph(AM);
IM = incidence(G);
IM = abs(IM);
AMline = round(IM'*IM) - 2*eye(numedges(G),numedges(G));
match = round(BK_MaxIS(AMline))>0;
ind = (sum(match,1) == k);
match = match(:,ind);
nTrans = size(match,2);

[s t] = find(bipart);
trans = zeros(k,nTrans);
for i = 1:nTrans
    trans(:,i) = s(match(:,i));
end
[];

