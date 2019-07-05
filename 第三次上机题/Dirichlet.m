n=20;
h = 1/n;
A = cell(n-1,n-1);   %系数矩阵
b = zeros((n-1)*(n-1),1);   %差分方程
B = zeros(n-1,n-1);
S = zeros(n-1,n-1);
%用于初始化S和B
for i = 1:n-1
        B(i,i) = -1/4;  % B = -1/4I
        S(i,i) = 1+h*h/4;   %S对角元均为1+h^2/4
        if (i>1)
            S(i,i-1) = -1/4;    %S次对角元均为-1/4
        end
        if (i<n-1)
            S(i,i+1) = -1/4;     %S次对角元均为-1/4
        end
end
%   初始化A和b
for i = 1:n-1
    for j = 1:n-1
        A{i,j} = zeros(n-1,n-1);
        b ((i-1)*(n-1)+j) = h*h/4*sin(i*h*j*h); 
        % b= h^2/4*f(i,j) + 1/4 (u(i+1,j)+u(i+1,j)+u(i+1,j)+u(i+1,j)) 
        %   b边界导入边界值
        if (i==1)
            b((i-1)*(n-1)+j) = b((i-1)*(n-1)+j) + j*h*j*h/4;  
        else if (i==n-1)
            b((i-1)*(n-1)+j) = b((i-1)*(n-1)+j) + j*h*j*h/4 + 1/4;   
            end
        end
        if (j==1)
            b((i-1)*(n-1)+j) = b((i-1)*(n-1)+j) + i*h*i*h/4;  
        else if (j==n-1)
            b((i-1)*(n-1)+j) = b((i-1)*(n-1)+j) + i*h*i*h/4 + 1/4;    
            end
        end
    end
    %替换A的主对角元和次对角元
    A{i,i} = S; %对角元为S
        %   次对角元为B
        if (i>1)
           A{i,i-1} = B;
        end
        if (i<n-1)
            A{i,i+1} = B;
        end
end
A = cell2mat(A);    %将A转换成普通矩阵

[x,r] = ConjugateGradient(zeros((n-1)*(n-1),1),A,b);
format short;
x=reshape(x,n-1,n-1)
format short e;
r=reshape(r,n-1,n-1)

X = zeros(n+1,n+1);
%   为x加入正方形区域的边界
for i = 0:n
    for j = 0:n
        if (i==0)
            X(i+1,j+1) = j*h*j*h;
            else if (i==n)
                    X(i+1,j+1) = j*h*j*h + 1;
                else if (j==0)
                        X(i+1,j+1) = i*h*i*h;
                     else if (j==n) 
                             X(i+1,j+1) = i*h*i*h + 1;
                         else
                                X(i+1,j+1) = x((i-1)*(n-1)+j);
                         end
                    end
                end
        end
    end
end
format short;
X
%surf(0:h:1,0:h:1,X);

