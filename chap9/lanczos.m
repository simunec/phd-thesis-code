function [V, T] = lanczos(A, b, m)
	% [V, T] = lanczos(A, b, m)
	% Runs the Lanczos algorithm for m iterations with matrix A and vector b 
	% Output:
	%	V 	n x m+1 matrix, orthonormal basis of K_{m+1}(A,b)  
	% 	T 	m+1 x m tridiagonal matrix, such that A * V(:,1:m) = V * T

n = length(b);
V = [b/norm(b), zeros(n,m)];
T = zeros(m+1,m);

w = A * V(:,1);
T(1,1) = w' * V(:,1);
w = w - T(1,1) * V(:,1);
T(2,1) = norm(w);
if m > 1
    T(1,2) = T(2,1);
end
V(:,2) = w / T(2,1);

for k = 2 : m
    w = A * V(:,k);
    T(k,k) = w' * V(:,k);
    w = w - T(k,k) * V(:,k) - T(k-1,k) * V(:,k-1);
    T(k+1,k) = norm(w);
    V(:,k+1) = w / T(k+1,k);
    
    if k < m
        T(k,k+1) = T(k+1,k);
    end
end
