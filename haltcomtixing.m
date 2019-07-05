function [T2n] = haltcomtixing( f,n,a,b,Tn )
%步长折半的复合梯形
%T2n = 1/2Tn + h/2*Σ(0,n-1)f(xk+1/2)
h = (b-a)/n;
T2n = 0;
for i = 0:n-1
    T2n = T2n+f(a+(i+1/2)*h);
end
T2n = 1/2*Tn + h/2*T2n;
end

