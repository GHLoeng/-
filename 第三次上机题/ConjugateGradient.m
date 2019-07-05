function [x,r] = ConjugateGradient(x,A,b )
%	�����ݶȷ�

r = b - A*x;
p = r;
k = 0;
while k==0 || max(abs(r))>=1e-4
    a = (r' * p) / (p' * A * p);    %������������
    x = x + a * p;   %���½�
    r = r - a * A * p;  %�����²в�����
    beta =  - (r' * A * p) / (p' * A * p);  
    p = r + beta * p;   %�����µ���������
    k = k+1;
end
k
end

