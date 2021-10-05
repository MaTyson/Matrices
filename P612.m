clc
clear

dim = [2 3 6 9 12 15 18 30 50 100 200 300 400 500];

for k = 1:numel(dim)
    n = dim(k);
    H = hilb(n);
    
    [Qc,Rc] = gsc(H);
    [Qm,Rm] = gsm(H);
    
    I = eye(n);
    nc = norm(I - Qc'*Qc,2);
    nm = norm(I - Qm'*Qm,2);
    
    plot(n,log10(nc),'*b');
    hold on
    plot(n,log10(nm),'dm');
    ylabel('log da norma','fontweight','bold','fontsize',16)
    xlabel('dim matriz','fontweight','bold','fontsize',16)
    title('QR matriz de Hilbert','fontweight','bold','fontsize',16)
    
    legend({'GSC', 'GSM'}, ...
    'Location', 'SouthEast','fontweight','bold','fontsize',12)
end