function x = Jacobi ( x,A,b )
%   �ſɱȵ�����
%   y = x
%   xi = (bi - ��(1,i-1)aijyj-��(i+1,n)aijyj)/aii

n = size(A,1);
y = x;
for i = 1:n
    sum = 0;
    for j = 1:i-1
        sum = A(i,j)*y(j) + sum;
    end
    for j = i+1:n
        sum = A(i,j)*y(j) + sum;
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
            sum = A(i,j)*y(j) + sum;
        end
        for j = i+1:n
            sum = A(i,j)*y(j) + sum;
        end
        x(i) = (b(i) - sum) / A(i,i)   ;
    end
    x
    InfiniteNorm(x-y)
end
end

