function x_inf = InfiniteNorm( x)
%	��x�������
x_inf = 0;
[a,b] = size(x);
for i = 1:a
    sum = 0;
    for j = 1:b
        sum = sum + abs(x(i,j));
    end
    if x_inf < sum
        x_inf = sum;
end
end

