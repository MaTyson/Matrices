function [A,b,alpha] = o3(A,b,i,k)
    [~,n] = size(A);
    candidato = abs(A(k,k));
    r = k;
    for l=k+1:n
        if( abs(A(l,k))>candidato)
            candidato = abs(A(l,k));
            r = l;
        end
    end
    if( r~=k)
        for j = k:n
            aux = A(k,j);
            A(k,j) = A(r,j);
            A(r,j) = aux;
        end
        aux = b(k);
        b(k) = b(r);
        b(r) = aux;
    end
    alpha = A(i,k)/A(k,k);
    A(i,k)=0;
    for j=k+1:n
        A(i,j) = A(i,j) - alpha*A(k,j);
    end
end