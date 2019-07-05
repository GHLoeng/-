function [ D ] = innerProduct( f,g)
D=0;
[n,m] = size(f);
for i = 1:m
    D = D + f(i)*g(i);
end
end

