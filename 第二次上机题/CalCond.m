function cond = CalCond( Hn)
% ����������
cond = OneNorm(Hn) * OneNorm(inv(Hn));
end

