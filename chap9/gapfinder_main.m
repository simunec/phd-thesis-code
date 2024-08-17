function [gaps, trest, bound_upper, bound_lower, gaps_rough, its_lanc] = gapfinder_main(A, mus, delta, theta, bound_type, d, c)
	% [gaps, trest, bound_upper, bound_lower, gaps_rough, its_lanc] = gapfinder_main(A, mus, delta, theta, bound_type, d, c)
	%
	% Finds gaps in the spectrum of A by computing trace(P_mu), for all mu in mus
	%
	% Input:
	%	A				matrix
	%	mus				vector of shifts for spectral projectors P_mu
	%	delta			failure probability
	%	theta			target relative gap width (if 0 < theta < 1)
	%		if theta > 1, it is considered as its_lanc, number of Lanczos iterations
	%	bound_type		type of bound for the Lanczos algorithm:
	%					"residue", a posteriori bound via residue theorem
	%					"diff", estimate via consecutive differences
	%					"diffsafe", safe version of "diff" estimate
	%	d				number of previous iterations used in Lanczos error
	%					bound (if bound_type = "diff" or "diffsafe"), or
	%					number of extra iterations for consecutive difference
	%					estimate (if bound_type = "residue") 
	%			[optional, default: 1]
	%	c				safety factor for consecutive difference estimate
	%			[optional, default: 1]
	%
	% Output:
	%	gaps			cell array with pairs [a_i, b_i] (extrema of found gaps)
	%	trest			trace estimates
	%	bound_upper		upper bounds on trace
	%	bound_lower		lower bounds on trace
	%	rough_gaps		optional, rough gaps computed with delta = 1/2

	if theta <= 1
		% use theta and delta to determine its_lanc 
		% gapfinder should detect all gaps with relative width larger than theta with probability at least 1-delta
		n = size(A, 1);
		C = 1 + (1 - theta)/(sqrt(pi*theta));
		epsilon = delta^2/exp(1);
		its_lanc = ceil(1 + log(2*C*n / epsilon) / log((1+theta)/(1-theta)));
		% fprintf("Lanczos iterations: %d\n", its_lanc);
	else
		its_lanc = ceil(theta);
	end

	if nargin < 7
		c = 1;
	end
	if nargin < 6
		d = 1;
	end
	if nargin < 5
		bound_type = "diffsafe";
	end

	N_hutch = 1;	% number of Hutchinson's sample vectors; theory recommends only 1 
	% % detect_gaps.m may not work as intended for N_hutch > 1

	[trest, bound_upper, bound_lower] = gapfinder_traceest(A, N_hutch, its_lanc, mus, bound_type, d, c);
	epsilon = delta^2/exp(1);
	[gaps, ~] = detect_gaps(mus, bound_upper, bound_lower, epsilon);

	% Optional, rough gaps:
	if nargout >= 5
		delta_rough = 0.5;
		epsilon_rough = delta_rough^2/exp(1);
		[gaps_rough, ~] = detect_gaps(mus, bound_upper, bound_lower, epsilon_rough);
	end

end