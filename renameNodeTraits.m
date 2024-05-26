function [AMtree_new,traits_new] = renameNodeTraits(AMtree,traits,nPat)
traits_new = traits;
traits_new(traits_new == 0) = (nPat+1):length(traits_new);
aux = [1:length(traits); traits_new]';
aux = sortrows(aux,2);
AMtree_new = AMtree(aux(:,1),aux(:,1));
traits_new = aux(:,2);
traits_new(traits_new > nPat) = 0;