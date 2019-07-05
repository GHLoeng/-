format long
A = [5 2 1;-1 4 2;2 -3 10];
b = [-12 20 3]';
x = [0 0 0]';
w = 0.9;
x = SOR(x,A,b,w);