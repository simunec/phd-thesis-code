function [qf, UU, vv] = matfun_proj_qf(A, b, f, V, T)
	% [qf, UU, vv] = matfun_proj_qf(A, b, f, V, T)
	% y = matfun_proj_qf(A, b, f, V, T)
	%* Computes an approximation to b'*f(A)*b by projecting on the Krylov subspace spanned by V
	% V is a matrix with orthonormal columns, and the approximation is given by 
	% qf = e_1'*f(T)*e_1, 		with T = V'*A*V
	% 	[we assume that V(:, 1) = b/norm(b)]
	% The function of the small matrix T is computed via eigendecomposition
	% f			can be a vector of functions
	%
	% y			final solution, i-th column is solution for the i-th function in f
	% UU, vv	orthogonal matrix and vector such that T = UU*diag(vv)*UU'

	n = size(A, 1);
	m = size(V, 2);
	
	% check if f is a function handle or a cell array
	if (isa(f, 'function_handle'))
		multifun = false;
		nfuns = 1;
	elseif iscell(f)
		multifun = true;
		nfuns = length(f);
	else
		error('function f is neither a function handle nor a cell array');
	end 

	% Compute final solution
	[UU, vv] = eig(T, 'vector');
	e1 = [1; zeros(m-1, 1)];
	if (~multifun)
		Fm = UU * diag(f(vv)) * UU';
		qf = (e1'*(Fm*e1)) * norm(b)^2;
	else
		qf = zeros(1, nfuns);
		c = UU' * (e1 * norm(b));
		for i = 1:nfuns
			qf(i) = c' * (f{i}(vv) .* c);
		end
	end

end
