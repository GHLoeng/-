format rat
A = [6  3;3 2];
b = [0 -1]';
x = [0 0]';
x = ConjugateGradient(x,A,b);