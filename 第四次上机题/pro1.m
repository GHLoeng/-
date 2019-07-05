x = [0 0.5 0.6 0.7 0.8 0.9 1.0];
y = [1 1.75 1.96 2.19 2.44 2.71 3.00];
y = y + 2*(rand(1,7)-0.5);
one = ones(1,7);
%y = a+bx
A = [innerProduct(one,one) innerProduct(one,x);
        innerProduct(x,one) innerProduct(x,x)];
b = [innerProduct(y,one);innerProduct(y,x)];
X = A\b;
a = X(1)
b = X(2)
s = 0;
for i = 1:7
    s = s + ((a+b*x(i))-y(i)).^2;
end
s
hold on
lx = 0:0.01:1;
plot(lx,a+b*lx);
scatter(x,y)






