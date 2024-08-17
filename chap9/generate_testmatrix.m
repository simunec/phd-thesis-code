function [A, U, D] = generate_testmatrix(ev)
	% function [A, U, D] = generate_testmatrix(ev)
	% Generates a tridiagonal test matrix A with eigenvalues ev 
	% A = U*D*U', 		D = diag(ev)

	n = length(ev);
	D = diag(ev);
	[Q, ~] = qr(randn(n));		% random orthogonal matrix
	M = Q*D*Q';
	[P, A] = hess(M);			% P*A*P' = M = Q*D*Q'
	A = tril(A,1);				% remove numerical zeroes
	A = 0.5*(A + A');			% ensure symmetry
	U = P'*Q;

end

