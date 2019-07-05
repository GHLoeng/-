function [x,r] = ConjugateGradient(x,A,b )
%	共轭梯度法

r = b - A*x;
p = r;
k = 0;
while k==0 || max(abs(r))>=1e-4
    a = (r' * p) / (p' * A * p);    %计算搜索步长
    x = x + a * p;   %更新解
    r = r - a * A * p;  %计算新残差向量
    beta =  - (r' * A * p) / (p' * A * p);  
    p = r + beta * p;   %计算新的搜索方向
    k = k+1;
end
k
end

