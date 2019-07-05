function [Sn,Rs] = comsimpson( f,N,a,b )
%   ¸´ºÏÐÁÆÕÉ­¹«Ê½
n = N/2;
h = (b-a)/n;
x = a:h:b;
Sn=f(a)+f(b);
for i = 0:n-1
    Sn = Sn + 4*f(a+h*i+h/2);
end
for i = 1:n-1
    Sn = Sn + 2*f(a+h*i);
end
Sn = Sn * h / 6;
syms x
f1 = f(fminbnd(matlabFunction(diff(sym(f),'x',4)),a,b));
f2 = f(fminbnd(matlabFunction(-diff(sym(f),'x',4)),a,b));
if f1>f2
    fmax = f1;
else
    fmax = f2;
end
Rs = h.^4/2880*(b-a)*fmax;
end

