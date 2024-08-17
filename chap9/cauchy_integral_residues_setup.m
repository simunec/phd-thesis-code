function [alpha, beta, gamma, theta, a2b2, abc2] = cauchy_integral_residues_setup(Xm, dm, hm)
	% [alpha, beta, gamma, theta, a2b2, abc2] = cauchy_integral_residues_setup(Xm, dm, hm)
	%* Performs some preliminary computations for function cauchy_integral_residues_quadform
	% we assume Lanczos method instead of rational Krylov, so Km = eye(m)
	% Am: m x m projected matrix (symmetric!)
	% [Xm, diag(dm)] = eig(Am);		% Am = Xm * diag(dm) * Xm^T
	% hm: last row of (m+1) x m  projected matrix Hm
	% Km: m x m projected matrix (assumed to be with m+1-th row = 0; in this case Km = eye(m))
	% Output: if Am = Xm diag(dm) Xm^{-1} and am^T = hm^T Km^{-1}
	% alpha = am^T * Xm
	% beta^T = Xm^{-1} * e1
	% gamma: coefficients obtained from alpha, beta and theta required for quadratic form bound
	% theta: eigenvalues of Am
	% a2b2 = alpha.^2 * beta.^2
	% abc2 = 2*alpha.*beta.*gamma
	

	m = length(dm);
	am = hm(:).';
	alpha = am * Xm;	alpha = alpha(:);
	beta = Xm \ [1; zeros(m-1, 1)];
	theta = dm;

	gamma = zeros(m, 1);
	for j = 1:m
		I = [1:j-1, j+1:m];
		gamma(j) = sum(alpha(I).*beta(I) ./ (theta(j) - theta(I)));
	end

	% Compute additional quantities:
	if nargout > 4
		ab = alpha.*beta;
		a2b2 = (ab).^2;
		abc2 = 2*ab.*gamma;
	end
end