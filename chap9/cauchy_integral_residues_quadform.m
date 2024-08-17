function [bound, lbound] = cauchy_integral_residues_quadform(a2b2, abc2, theta, f, df, lambda)
	% [bound, lbound] = cauchy_integral_residues_quadform(a2b2, abc2, theta, f, df, lambda)
	%* Computes bound for approximation of b'*f(A)*b after m iterations of rational Krylov method with projected matrix Am
	% uses a more efficient implementation for several functions f
	% requires calling cauchy_integral_residues_setup with all output arguments
	%
	% Input:
	%	a2b2: alpha.^2 * beta.^2
	%	abc2: 2*alpha.*beta.*gamma
	%	theta: eigenvalues of projected matrix Am
	%		[a2b2, abc2, theta must all be column vectors]
	%	f: functions f (cell array)
	%	df: derivative of f (cell array)
	%		[f and df should allow vector input]
	%	lambda: candidate eigenvalues, e.g. discretization of spectral interval [lmin, lmax]
	%		[row vector]
	%
	% Output:
	%	bound: upper bound for the error
	%	lbound: lower bound (optional, usually less accurate)	
	%
	% Can be slow if length(lambda) and length(f) are large
	% There may be issues when Am has some repeated or very close eigenvalues 
	% (NaN and Inf appear already in the gamma coeffs computed in setup phase)

	lambda = lambda(:).';
	N = length(lambda);
	m = length(theta);
	nfuns = length(f);
	
    % Make sure that we deal with column vectors
    a2b2 = a2b2(:);
    abc2 = abc2(:);
    theta = theta(:);
	
	% Precompute function values:
	flambda = zeros(1, N, nfuns);
	ftheta = zeros(m, 1, nfuns);
	dftheta = zeros(m, 1, nfuns);
	for j = 1:nfuns
		flambda(1, :, j) = f{j}(lambda);
		ftheta(:, 1, j) = f{j}(theta);
		dftheta(:, 1, j) = df{j}(theta);
	end
	inv_lambda_m_theta = 1 ./ (lambda - theta);
	%! Trick to efficiently ignore NaNs in the sum below:
	%! only works for functions with numerator equal to zero when denominator is infinite (such as sign function)
	inv_lambda_m_theta(isinf(inv_lambda_m_theta)) = 1;	
	inv_lambda_m_theta2 = inv_lambda_m_theta.^2;

	bound = zeros(1, nfuns);
	if nargout > 1
		lbound = zeros(1, nfuns);
	end
	for j = 1:nfuns
        gm_vec1 = a2b2 .* ( (flambda(1, :, j) - ftheta(:, 1, j)).*inv_lambda_m_theta2 - dftheta(:, 1, j).*inv_lambda_m_theta);
		gm_vec2 = abc2 .*  (flambda(1, :, j) - ftheta(:, 1, j)).*inv_lambda_m_theta;
		gm = sum(gm_vec1, 1) + sum(gm_vec2, 1); 
		% gm = sum(gm_vec1, 1, 'omitnan') + sum(gm_vec2, 1, 'omitnan');  % not required if using trick above

		bound(j) = max(abs(gm));
		if nargout > 1
			lbound(j) = min(abs(gm));
			if (min(gm) < 0 && max(gm) > 0)
				lbound(j) = 0;
			end
		end
	
	end

end