function [x,r] = UsefulConjugateGradient(x,A,b )
%	ʵ�ù����ݶȷ�

r = b - A*x;
p = r;
k = 0;
while k == 0 || max(abs(r))>=1e-4
    a = (r' * p) / (p' * A * p);    %������������
    x = x + a * p;   %���½�
    rp = r; %������һ���в�����
    r = r - a * A * p;  %���²в�����
    beta =  - (r' * r) / (rp' * rp);  
    p = r + beta * p;   %�����µ���������
    k = k+1;
end
k
end

