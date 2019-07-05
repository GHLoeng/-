function cutoff = standardChop(coeffs, tol)
%STANDARDCHOP A sequence chopping rule of ''standard'' (as opposed to ''loose'' or
% ''strict'') type, that is, with an input tolerance TOL that is applied with some
% flexibility. This code is used in all parts of Chebfun that make chopping
% decisions, including chebfun construction (CHEBTECH, TRIGTECH), solution of
% ODE BVPs (SOLVEBVP), solution of ODE IVPs (ODESOL), simplification of chebfuns % (SIMPLIFY), and Chebfun2. See J. L. Aurentz and L. N. Trefethen, ''Chopping a % Chebyshev series'', arXiv, December 2015.
%
% Input:
%
% COEFFS实数或复数的非空行列向量，通常为切比雪夫或傅立叶系数。
%
% TOL（0,1）中代表目标相对准确度的数字。
%	TOL通常将设置为Chebfun CHEBFUNEPS参数，有时乘以诸如vglobal / vlocal之类的因子来构建局部全局chebfuns。
%	默认值：机器epsilon（MATLAB EPS）。
%
% Output:
%
% CUTOFF A positive integer.
%	If CUTOFF == length(COEFFS), then we are ''not happy'':
%	a satisfactory chopping point has not been found.
%	If CUTOFF < length(COEFFS), we are ''happy'' and CUTOFF
%	represents the last index of COEFFS that should be retained.
%
% Examples:
%
% coeffs = 10.^-(1:50); random = cos((1:50).^2);
% standardChop(coeffs) % = 18
% standardChop(coeffs + 1e-16*random) % = 15
% standardChop(coeffs + 1e-13*random) % = 13
% standardChop(coeffs + 1e-10*random) % = 50
% standardChop(coeffs + 1e-10*random, 1e-10) % = 10
% Jared Aurentz and Nick Trefethen, July 2015.
%
% Copyright 2015 by The University of Oxford and The Chebfun Developers.  See http://www.chebfun.org/ for Chebfun information.
% STANDARDCHOP normally chops COEFFS at a point beyond which it is smaller than
% TOL^(2/3). COEFFS will never be chopped unless it is of length at least 17 and
% falls at least below TOL^(1/3). It will always be chopped if it has a long
% enough final segment below TOL, and the final entry COEFFS(CUTOFF) will never
% be smaller than TOL^(7/6). All these statements are relative to
% MAX(ABS(COEFFS)) and assume CUTOFF > 1. These parameters result from
% extensive experimentation involving functions such as those presented in  the paper cited above. They are not derived from first principles and  there is no claim that they are optimal.
% Set default if fewer than 2 inputs are supplied:

%STANDARDCHOP通常会在比TOL ^（2/3）更小的点上截断COEFFS。
%除非长度至少为17并且至少低于TOL ^（1/3），否则COEFFS将永远不会被截断。
%如果它在TOL下有足够长的最后段，并且最终的元素COEFFS（CUTOFF）永远不会小于TOL ^（7/6），它将总是被切断。
%所有这些陈述都是相对的
%MAX（ABS（COEFFS））并且假定CUTOFF> 1。
%这些参数由涉及诸如上文引用的论文中的那些功能的广泛实验产生。
%它们不是源于第一原则，也没有声称它们是最优的。
%如果提供的输入少于2个，则设置为默认值：

if ( nargin < 2 )
    p = chebfunpref; 
    tol = p.chebfuneps;
end
% 检查TOL的大小：
if ( tol >= 1 ) 
    cutoff = 1;
    return
end
% 确保COEFFS的长度至少为17：
    n = length(coeffs); 
    cutoff = n;
if ( n < 17 ) 
    return
end
%步骤1：将COEFFS转换为一个新的单调不增加的向量ENVELOPE，并将其归一化为以值1开头。
b = abs(coeffs); 
m = b(end)*ones(n, 1); 
for j = n-1:-1:1 
    m(j) = max(b(j), m(j+1));
end
if ( m(1) == 0 ) 
    cutoff = 1; 
    return
end
envelope = m/m(1);
% Step 2: Scan ENVELOPE for a value PLATEAUPOINT, the first point J-1, if any,
% that is followed by a plateau. A plateau is a stretch of coefficients
% ENVELOPE(J),...,ENVELOPE(J2), J2 = round(1.25*J+5) <= N, with the property
% that ENVELOPE(J2)/ENVELOPE(J) > R. The number R ranges from R = 0 if
% ENVELOPE(J) = TOL up to R = 1 if ENVELOPE(J) = TOL^(2/3). Thus a potential
% plateau whose starting value is ENVELOPE(J) ~ TOL^(2/3) has to be perfectly
% flat to count, whereas with ENVELOPE(J) ~ TOL it doesn't have to be flat at
% all. If a plateau point is found, then we know we are going to chop the % vector, but the precise chopping point CUTOFF still remains to be determined % in Step 3.

%步骤2：扫描ENVELOPE得到一个值PLATEAUPOINT，第一个点J-1（如果有的话），后面跟着平稳水平。 
%平稳水平是一系列ENVELOPE(J),...,ENVELOPE(J2), J2 = round(1.25*J+5) <= N,具有ENVELOPE(J2)/ENVELOPE(J) > R
%如果ENVELOPE（J）= TOL，R = 0，如果ENVELOPE（J）= TOL ^（2/3），R = 1。 
%因此，一个初始值为ENVELOPE（J）?TOL ^（2/3）的潜在平稳水平必须完全平坦才能计数，
%而使用ENVELOPE（J）?TOL时，它不一定是平坦的。
%如果找到一个平稳点，那么我们知道我们要截断矢量，但精确的截断点CUTOFF仍然在步骤3中确定。
for j = 2:n
    j2 = round(1.25*j + 5);
    if ( j2 > n )
    % there is no plateau: exit 
    return
    end
    e1 = envelope(j); 
    e2 = envelope(j2); 
    r = 3*(1 - log(e1)/log(tol)); 
    plateau = (e1 == 0) | (e2/e1 > r);
    if ( plateau )
        % a plateau has been found: go to Step 3
        plateauPoint = j - 1;
        break
    end
end
% Step 3: fix CUTOFF at a point where ENVELOPE, plus a linear function % included to bias the result towards the left end, is minimal.
%
% Some explanation is needed here. One might imagine that if a plateau is
% found, then one should simply set CUTOFF = PLATEAUPOINT and be done, without % the need for a Step 3. However, sometimes CUTOFF should be smaller or larger % than PLATEAUPOINT, and that is what Step 3 achieves.
%
% CUTOFF should be smaller than PLATEAUPOINT if the last few coefficients made
% negligible improvement but just managed to bring the vector ENVELOPE below the
% level TOL^(2/3), above which no plateau will ever be detected. This part of % the code is important for avoiding situations where a coefficient vector is % chopped at a point that looks ''obviously wrong'' with PLOTCOEFFS.
%
% CUTOFF should be larger than PLATEAUPOINT if, although a plateau has been
% found, one can nevertheless reduce the amplitude of the coefficients a good % deal further by taking more of them. This will happen most often when a
% plateau is detected at an amplitude close to TOL, because in this case, the
% ''plateau'' need not be very flat. This part of the code is important to % getting an extra digit or two beyond the minimal prescribed accuracy when it % is easy to do so.

%第3步：修正CUTOFF在一个点，其中ENVELOPE，加上一个线性函数，将结果偏向左端，是最小的。
%这里需要一些解释。可以想象，如果找到平稳水平，那么应该简单地设置CUTOFF = PLATEAUPOINT并完成，而不需要第3步。
%但是，有时CUTOFF应该小于或大于PLATEAUPOINT，这就是步骤3实现的目标。
%?如果最后几个系数的改进可以忽略不计，但只是设法使矢量ENVELOPE低于TOL ^（2/3）的水平，CUTOFF应小于PLATEAUPOINT
% 高于此水平将不会检测到平稳水平。这部分代码对于避免系数向量在PLOTCOEFFS看起来“明显错误”的点处被截断的情况很重要。
%?如果虽然已经找到了平稳水平，但仍然可以通过采用更多的系数来进一步降低系数的幅度，所以CUTOFF应该大于PLATEAUPOINT。
%这通常会发生在接近TOL的幅度处检测到平稳水平，因为在这种情况下，“平稳水平”不需要非常平坦。
%这部分代码对于获得超出最小规定精度的额外数字十分重要，因为这很容易实现。
if ( envelope(plateauPoint) == 0 ) 
    cutoff = plateauPoint;
else
    j3 = sum(envelope >= tol^(7/6)); 
    if ( j3 < j2 ) 
        j2 = j3 + 1;
        envelope(j2) = tol^(7/6);
    end
    cc = log10(envelope(1:j2)); 
    cc = cc(:);
    cc = cc + linspace(0, (-1/3)*log10(tol), j2)';
    [~, d] = min(cc); 
    cutoff = max(d - 1, 1); 
end

end
