function [Q,R] = gsm(A)
    [m,n] = size(A);
    Q = zeros(m,n);
    R = zeros(n,n);
    v = zeros(m,n);
    
    for j = 1:n
        v(:,j) = A(:,j);
    end
    for k = 1:n
        R(k,k) = norm(v(:,k),2);
        Q(:,k) = v(:,k)/R(k,k);
        for j = k+1:n
            R(k,j) = Q(:,k)'*v(:,j);
            v(:,j) = v(:,j) - R(k,j)*Q(:,k);
        end
    end
end
    
