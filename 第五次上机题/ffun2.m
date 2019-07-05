function [ dy ] = ffun2( x,y )
v1 = 1;
v2 = 2;
dy = zeros(2,1);
dy(1) = y(2);
dy(2) = v1/v2/x*sqrt(1+y(2).^2);
end

