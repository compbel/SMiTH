function [AM,patients,seqID] = phytree2graph(Z,delimeter,pos)

names = get(Z,'LeafNames');
nSeq = length(names);
patients = zeros(1,2*nSeq-1);
for i = 1:nSeq
    C = strsplit(names{i},delimeter);
    patients(i) = str2num(C{pos})+1;
end

[A,~,~] =  getmatrix(Z);
AM = full(A);

