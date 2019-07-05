function x =  SteepestDescent( x ,A, b)
%	最速下降法
%   r(k+1) = b - Ax(k+1)
%   x(k+1) = xk + ark

r = b - A * x;
k = 0;
while (k<4)
    a = (r' * r) / ( r' * A * r);
    x = x + a * r
    r = r - a * A * r;
    k = k+1;
end

end

