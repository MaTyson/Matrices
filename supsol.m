function x = supsol(U,b)
    [~,n] = size(U);
    x = zeros(n,1);
    for k = n:-1:1
        s = 0;
        for j = k+1:n
            s = s + U(k,j)*x(j);
        end
        x(k) = (b(k) - s)/U(k,k);
    end
end