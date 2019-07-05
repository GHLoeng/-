function [Tn,Rt] = comtixing( f,N,a,b )
    %   复合梯形公式
n = N;
h = (b-a)/n;
x = a:h:b;
Tn=f(a)+f(b);
for i = 1:n-1
    Tn = Tn + 2*f(a+h*i);
end
Tn = Tn * h / 2;
syms x
f1 = f(fminbnd(matlabFunction(diff(sym(f),'x',2)),a,b));
f2 = f(fminbnd(matlabFunction(-diff(sym(f),'x',2)),a,b));
if f1>f2
    fmax = f1;
else
    fmax = f2;
end
Rt = h.^2/2*(b-a)*fmax;
end

