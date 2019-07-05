function x = G_S( x,A,b )
%   G-Sµü´ú·¨
%   xi = (bi - ¦²(1,i-1)aijxj-¦²(i+1,n)aijxj)/aii

n = size(A,1);
y = x;
for i = 1:n
    sum = 0;
    for j = 1:i-1
        sum = A(i,j)*x(j) + sum;
    end
    for j = i+1:n
        sum = A(i,j)*x(j) + sum;
    end
    x(i) = (b(i) - sum) / A(i,i) ;
end
x
InfiniteNorm(x-y)
while (InfiniteNorm(x-y) >= 1e-2)
    y = x;
    for i = 1:n
        sum = 0;
        for j = 1:i-1
            sum = A(i,j)*x(j) + sum;
        end
        for j = i+1:n
            sum = A(i,j)*x(j) + sum;
        end
        x(i) = (b(i) - sum) / A(i,i)   ;
    end
    x
    InfiniteNorm(x-y)
end

end

