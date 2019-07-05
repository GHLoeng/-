n = 20;
hold on
for n = n:20
t = [-1:2/n:1-2/n]';
f = exp(t);
for l = 3:2:9
    A = zeros( n,l+1);
    for i = 1:n
        for j = 1:l+1
            A(i,j) = legendre(j-1,t(i));
        end
    end
    At = A';
    G = At*A;
    b = At*f;
    L = Cholesky(G,l+1);
    y = L\b;
    x = L'\y;
    q = A*x-f;
    s = 0;
    for i = 1:n
        s = s + q(i)*q(i);
    end
    s
    xlabel('n');
    ylabel('ln(cond(A))');
    scatter(n,log(cond(A)));
end
end