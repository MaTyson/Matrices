function [Q,R] = gsc(A)
    [m,n] = size(A);
    Q = zeros(m,n);
    R = zeros(n,n);
    
    for k = 1:n
        v = A(:,k);
        for j = 1:k-1
            R(j,k) = Q(:,j)'*A(:,k);
            v = v - R(j,k)*Q(:,j);
        end
        R(k,k) = norm(v,2);
        Q(:,k) = v/R(k,k);
    end
end