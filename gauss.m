function [A,b] = gauss(A,b)
    [~,n] = size(A);
    for k = 1:n-1
        for i = (k+1):n
            [A,b,alpha] = o3(A,b,i,k);
            b(i) = b(i) -alpha*b(k);
        end
    end
end

