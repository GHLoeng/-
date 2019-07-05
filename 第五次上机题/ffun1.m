function [ dy ] = ffun1( t,y )
v1 = 1;
v2 = 2;
dy = zeros(2,1);
dy(1) = v2/(sqrt(1+((y(2)-v1*t)/y(1)).^2));
dy(2) = v2/(sqrt(1+((y(2)-v1*t)/y(1)).^2)) * ((y(2)-v1*t)/y(1)).^2;
end

