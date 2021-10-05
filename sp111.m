clc
clear

n = 3000;
ds = 0.01;
A = 20*sprand(n,n,ds,0.1);
A = A - 10*spones(A);
xsol = 4*rand(n,1) - 2;
b = A*xsol;

[L,U,P] = lu(A);
c = P*b;
y = L\c;
x = U\y;
norm(x-xsol,2)

p2 = colperm(A);
B = A(:,p2);
[LB,UB,PB] = lu(B);
c = PB*b;
y = LB\c;
z = UB\y;
x(p2) = z;
norm(x-xsol,2)
