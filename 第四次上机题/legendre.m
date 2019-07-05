function [ p ] = legendre ( n,t )
%   勒让德多项式
%   n阶数，t自变量
if n == 0
    p = 1;
else if n == 1
        p = t;
    else
        p = ((2*n-1)*t*legendre(n-1,t)-(n-1)*legendre(n-2,t))/n;
    end
end
end

