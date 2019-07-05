function cond = CalCond( Hn)
% 计算条件数
cond = OneNorm(Hn) * OneNorm(inv(Hn));
end

