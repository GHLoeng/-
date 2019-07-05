n = input ('����������ά��n= ');
Hn = BuildHilbert(n);
x_exact = ones(n,1);
b = Hn * x_exact;
format rat
L = Cholesky(Hn,n);
LT = L.';
y = L\b; % Ly=b
x_appro = LT\y ;% LTx=y


r = b - Hn * x_appro; % r = b - Hnx'
sx = x_appro - x_exact; % ��x = x' - x

ri = InfiniteNorm(r);
sxi = InfiniteNorm(sx);
cond = CalCond(Hn)

fprintf('r���������%e\n',ri);
fprintf('��x���������%e\n',sxi);