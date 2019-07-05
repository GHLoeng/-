function x = ConjugateGradient(x,A,b )
%	¹²éîÌİ¶È·¨

r = b - A*x;
p = r;
k = 0;
while k<4
    a = (r' * p) / (p' * A * p);
    x = x + a * p
    r = r - a * A * p;
    beta =  - (r' * A * p) / (p' * A * p);
    p = r + beta * p;
    k = k+1;
end
end

