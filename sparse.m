clear
clc

dim = [1000 2000 3000 4000 5000 6000 7000 8000 9000 10000];
tempo = zeros(length(dim),3);
for k = 1:length(dim)
    t0 = tic;
    n = dim(k)
    A = 20*sprand(n,n,0.001,1);
    A = A - 10*spones(A);
    
    b = 10*rand(n,1) - 5;
    
    t1 = tic;
    x = A\b;
    tempo(k,1) = toc(t1);
    
    C = full(A);
    
    t2 = tic;
    w = C\b;
    tempo(k,2) = toc(t2);
    
    tempo(k,3) = toc(t0);
    total = tempo(k,3)
end
tempo