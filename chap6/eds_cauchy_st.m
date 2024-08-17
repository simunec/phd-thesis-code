function [poles, bound] = eds_cauchy_st(a, b, m, z)
	%* Poles for (m, m) rational approximation on [a, b] of Cauchy-Stieltjes functions
	% Asymptotically optimal poles in (-Inf, 0) constructed with an EDS
	% Algorithm taken from [1]
	% z can be a custom (irrational) starting value for the EDS in [0, 1]
	% bound contains the bound (without constant) taken from [1, Table 1]
	% Cauchy-Stieltjes functions (a subclass of Laplace-Stieltjes) can be written as integral:
	%	f(z) = \int_0^\infty {t + z)^{-1} \mu(t) dt, 	z > 0
	% Now fixed by removing cancellation for large values of b/a

	% [1] Massei, Robol, "Rational Krylov for Stieltjes 
	% matrix functions: convergence and pole selection", BIT, 2020

	%* Magic to transform from [a, b] U (-Inf, 0] to [-1, -ahat] U [ahat, 1]:
	Delta = sqrt(b^2 - a*b);
	ahat_old = (b - Delta)/(b + Delta);
	ahat = a / (sqrt(b) + sqrt(b-a))^2;	% to avoid cancellation
	
	%* Construct equidistributed sequence s on [0, 1]:
	if nargin < 4
		z = 1/sqrt(2);	
	end
	s = zeros(1, m);

	for k = 1:m
		s(k) = k*z - floor(k*z);
	end

	%* Evaluate elliptic function on s:
	[~, Kp] = ellipk(ahat);
	[~, ~, dn] = ellipj(s * Kp, 1 - ahat^2);
	poles = -dn;		% poles on [-1, -ahat]
	poles_old = ((b + Delta)*poles + (b - Delta)) ./ (1+poles);		% poles on (-Inf, 0]
	poles = ((b + Delta)*poles +  sqrt(b)*a / (sqrt(b) + sqrt(b-a))) ./ (1+poles);		% poles on (-Inf, 0], avoiding cancellation

	%* Convergence rate:
	R = exp(-pi^2 / log(4*4*b/a));
	bound = R.^[1:m];

	% %* Poles on [-b, -a] for Laplace-Stieltjes functions:
	% %* For Laplace-Stieltjes we should have poles on [-b, -a]:
	% poles_lap = -dn;		% poles on [-1, -ahat]
	% poles_lap = poles_lap * (b-a)/(1-ahat) + (a - ahat*(b-a)/(1-ahat)) ;

	%* Zolotarev optimal poles should be:
	% zz = ((2*m + 1) - [1:m]) / (2*m);
	% zolpoles = ellipj(zz*Kp, 1 - ahat^2);

end