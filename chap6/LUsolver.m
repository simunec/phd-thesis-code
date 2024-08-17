function y = LUsolver(L, U, w, perm)
% y = LUsolver(L,U,w, perm) 
% Computes the solution to the linear system A y = w, 
% using the LU factorization A(perm, perm) = L*U
% the solution is computed as y(perm) = U \ (L \ w(perm))
% using back- and forward-substitution
% U is upper triangular
% L is lower triangular

	perm = perm(:);
	y(perm) = U \ (L \ w(perm));	
	y = y(:);

end
