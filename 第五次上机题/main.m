[x,y]=ode45(@ffun2,[100 1e-15],[0 0]);
plot(x,y(:,1))
xlabel('x/m')
ylabel('y/m')
a = 100;
v1 = 1;
v2 = 2;
y(end/2)    %����ֵ
yp = a*v1*v2/(v2*v2-v1*v1)  %��ȷֵ
e = y(end/2) - yp   %���