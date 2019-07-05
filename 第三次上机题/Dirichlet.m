n=20;
h = 1/n;
A = cell(n-1,n-1);   %ϵ������
b = zeros((n-1)*(n-1),1);   %��ַ���
B = zeros(n-1,n-1);
S = zeros(n-1,n-1);
%���ڳ�ʼ��S��B
for i = 1:n-1
        B(i,i) = -1/4;  % B = -1/4I
        S(i,i) = 1+h*h/4;   %S�Խ�Ԫ��Ϊ1+h^2/4
        if (i>1)
            S(i,i-1) = -1/4;    %S�ζԽ�Ԫ��Ϊ-1/4
        end
        if (i<n-1)
            S(i,i+1) = -1/4;     %S�ζԽ�Ԫ��Ϊ-1/4
        end
end
%   ��ʼ��A��b
for i = 1:n-1
    for j = 1:n-1
        A{i,j} = zeros(n-1,n-1);
        b ((i-1)*(n-1)+j) = h*h/4*sin(i*h*j*h); 
        % b= h^2/4*f(i,j) + 1/4 (u(i+1,j)+u(i+1,j)+u(i+1,j)+u(i+1,j)) 
        %   b�߽絼��߽�ֵ
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
    %�滻A�����Խ�Ԫ�ʹζԽ�Ԫ
    A{i,i} = S; %�Խ�ԪΪS
        %   �ζԽ�ԪΪB
        if (i>1)
           A{i,i-1} = B;
        end
        if (i<n-1)
            A{i,i+1} = B;
        end
end
A = cell2mat(A);    %��Aת������ͨ����

[x,r] = ConjugateGradient(zeros((n-1)*(n-1),1),A,b);
format short;
x=reshape(x,n-1,n-1)
format short e;
r=reshape(r,n-1,n-1)

X = zeros(n+1,n+1);
%   Ϊx��������������ı߽�
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

