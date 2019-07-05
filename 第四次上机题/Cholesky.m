function A = Cholesky( A , n )
%   Cholesky 分解算法
%   输入：A，n 输出：A
%   (1) L11 = sqrt(a11)
%   (2) Li1(i>=2) = ai1/L11
%   (3) Ljj = sqrt(ajj - Σ(1,j-1)Ljk^2
%   (4) Lij(i>j) = (aij - Σ(1,j-1)LikLjk) / Ljj

for j = 1:n
    for k = 1:j-1
        A(j,j) = A(j,j) - A(j,k)*A(j,k);
    end
    
    A(j,j) = sqrt(A(j,j)); %(3)
    
    for i = j+1:n
        for k = 1:j-1
            A(i,j) = A(i,j) - A(i,k)*A(j,k);
        end
        A(i,j) = A(i,j) / A(j,j); %(4)
    end
end
for i = 1:n
    for j = i+1:n
        A(i,j) = 0;
    end
end
end

