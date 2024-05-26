function G = randTreeUniform(n)

    prufer = randi(n,1,n-2);
    deg = ones(1,n);
    for j = prufer
        deg(j) = deg(j)+1;
    end
    AM_G = zeros(n,n);
    for i = prufer
        j = find(deg==1,1);
        AM_G(i,j)=1;
        AM_G(j,i)=1;
        deg(i) = deg(i)-1;
        deg(j) = deg(j)-1;
    end

    pair = find(deg==1,2);
    AM_G(pair(1),pair(2))=1;
    AM_G(pair(2),pair(1))=1;

    G = graph(AM_G);
