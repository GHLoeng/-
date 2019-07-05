function [ R2n ] = romberg( f,a,b,e )
    %龙贝格求积公式
R = zeros(2,2);
n=2;
h = b-a;
R(1,1) = h/2*(f(a)+f(b));
while 1
    R(1,n) = haltcomtixing(f,2.^(n-2),a,b,R(1,n-1));
    for j = 2:n
        R(j,n-j+1) = ((4^(j-1))*R(j-1,n-j+2)-R(j-1,n-j+1))/((4^(j-1))-1);
    end
    if abs(R(n,1)-R(n-1,1))<e
        break;
    end
    n = n+1;
end
R2n = R(n,1);
R
end

