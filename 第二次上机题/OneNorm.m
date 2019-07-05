function x_1 = OneNorm( x )
%	ÇóxµÄ1-·¶Êý
x_1 = 0;
[a,b] = size(x);
for j = 1:a
    sum = 0;
    for i = 1:b
        sum = sum + abs(x(i,j));
    end
    if x_1 < sum
        x_1 = sum;
end
end

