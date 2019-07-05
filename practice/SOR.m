function x = SOR( x,A,b,w )
%   SORµü´ú·¨
%   xi = (1-w)xi + w(bi - ¦²(1,i-1)aijxj-¦²(i+1,n)aijxj)/aii

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
        x(i) = (1-w)*x(i) + w * (b(i) - sum) / A(i,i)   ;
    end
    x
    InfiniteNorm(x-y)
end
end

