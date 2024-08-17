function y = fracexpmvec(A, alpha, t, b)
% y = fracexpmvec(A, alpha, t, b) returns y = exp(-t A^alpha) b

	% The following algorithm is only stable for symmetric or well conditioned A
	[Q,D] = eig(A);
	y = Q* (diag( exp(-t*diag(D).^(alpha)) )*(Q\b));
	y = real(y);

end