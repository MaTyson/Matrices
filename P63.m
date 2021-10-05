clear
clc

%A = [5 2 1; 8 5 0; -1 1 -1];
%b = [3;7;2];
%xsol = A\b;

n = input('manda  dim ---> ')
A = round(10*rand(n,n) + 5);
%A = hilb(n);
xsol = round(2*rand(n,1)+1);
b = A*xsol;

ncond = cond(A,inf);
%[L, U, P] = lu(A);
[L, U, p] = lu(A,'vector');
%z = P*b;
z = b(p);
y = L\z;
x = U\y;
r = b - A*x;
nr = norm(r,inf)/norm(b,inf);
ex = xsol - x;
ne = norm(ex,inf)/norm(xsol,inf);
t = ne/norm(x,inf);
q = ncond*nr;
ncond
nr
ne
t
q
t<=q