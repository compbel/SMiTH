function divers = getDiversity(seqFile,separator,token)
seqF = fastaread(seqFile);
nSeq = length(seqF);
nm = char(seqF.Header);
seq = char(seqF.Sequence);
seqPatID = zeros(1,nSeq);
L = size(seq,2);
for i = 1:nSeq
    tokens = strsplit(nm(i,:), separator);
    seqPatID(i) = str2num(tokens{token})+1;
end
patientList = sort(unique(seqPatID));
nPat = length(patientList);


maxDist = zeros(1,nPat);
meanDist = zeros(1,nPat);
meanEntropy = zeros(1,nPat);
for i = 1:nPat
    ind = (seqPatID == patientList(i));
    seq_i = seq(ind,:);
    
    if sum(ind) > 1
        DM = pdist(seq_i, 'hamming');
        maxDist(i) = max(DM(:));
        meanDist(i) = mean(DM(:));
    else
        maxDist(i) = 0;
        meanDist(i) = 0;
    end

    totalEntropy = 0;
    for col = 1:L
        [~, ~, idx] = unique(seq_i(:, col));
        charCounts = accumarray(idx, 1);
        probabilities = charCounts / sum(charCounts);
        entropy = -sum(probabilities .* log2(probabilities));
        totalEntropy = totalEntropy + entropy;
    end
    meanEntropy(i) = totalEntropy / L;
end
divers = meanEntropy;
