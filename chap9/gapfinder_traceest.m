function [trest, bound_upper, bound_lower] = gapfinder_traceest(A, N_hutch, its_lanc, mus, bound_type, d, c)
	% [trest, bound_upper, bound_lower] = gapfinder_traceest(A, N_hutch, its_lanc, mus, bound_type, d, c)
	% Estimates trace(P_mu) for finding gaps in the spectrum of A, for all mu in mus

	% Input:
	%	A				matrix
	%	N_hutch			number of sample vectors for Hutchinson's algorithm 
	%			[recommended: 1]
	%	its_lanc		number of lanczos iterations
	%	mus				vector of shifts for spectral projectors P_mu
	%	bound_type		type of bound for the Lanczos algorithm:
	%					"residue", a posteriori bound via residue theorem
	%					"diff", estimate via consecutive differences
	%					"diffsafe", safe version of "diff" estimate
	%			[optional, default: "diffsafe"]
	%	d				number of previous iterations used in Lanczos error
	%					bound (if bound_type = "diff" or "diffsafe"), or
	%					number of extra iterations for consecutive difference
	%					estimate (if bound_type = "residue") 
	%			[optional, default: 1]
	%	c				safety factor for consecutive difference estimate
	%			[optional, default: 1]

	% Output:
	%	trest			trace estimates
	%	bound_upper		upper bounds on trace
	%	bound_lower		lower bounds on trace

	
	if nargin < 7
		c = 1;
	end
	if nargin < 6
		d = 1;
	end
	if nargin < 5
		bound_type = "diffsafe";
	end

	n = size(A, 1);
	nfuns = length(mus);
	fmu = cell(nfuns, 1);
	dfmu = cell(nfuns, 1);
	for k = 1:nfuns
		fmu{k} = @(z) (real(z) <= mus(k));		% step function
		dfmu{k} = @(z) 0;
	end		

	Omega = randn(n, N_hutch);					% Gaussian
	trest = zeros(1, nfuns);
	bound_upper = zeros(1, nfuns);
	bound_lower = zeros(1, nfuns);
	for k = 1:N_hutch
		if bound_type == "residue"
			[V, T] = lanczos(A, Omega(:, k), its_lanc);
		else
			[V, T] = lanczos(A, Omega(:, k), its_lanc+1);
		end

		if bound_type == "residue"
			% Compute trace approximation bound with residue theorem:
			bound_upper_cur = Inf*ones(1, nfuns);
			bound_lower_cur = zeros(1, nfuns);		
			for il = its_lanc-d+1:its_lanc
				% Compute trace approximation:
				[trestk, UU, vv] = matfun_proj_qf(A, Omega(:, k), fmu, V(:, 1:il), T(1:il, 1:il));	
				% Compute trace approximation bound:
				[~, ~, ~, theta, a2b2, abc2] = cauchy_integral_residues_setup(UU, vv, T(il+1, 1:il));
				lambda = linspace(min(theta), max(theta), 1000);
				[boundk, ~] = cauchy_integral_residues_quadform(a2b2, abc2, theta, fmu, dfmu, lambda);
				boundk = boundk * norm(Omega(:, k))^2;		% multiply by norm(Omega(:, k)^2 ~ n
	
				% Compute upper and lower bounds for Omega(:,k)' * P_mu * Omega(:,k):
				boundk_upper = cummin(trestk + boundk, 'reverse');
				boundk_lower = cummax(trestk - boundk);
	
				% Update bounds from last few iterations:
				bound_upper_cur = min(bound_upper_cur, boundk_upper);
				bound_lower_cur = max(bound_lower_cur, boundk_lower);
			end


		else
			% Compute trace approximation bound with consecutive difference:
			trestk = zeros(d+1, nfuns);
			if bound_type == "diff"
				bound_lower_cur = zeros(1, nfuns);		
				bound_upper_cur = Inf*ones(1, nfuns);
			elseif bound_type == "diffsafe"
				bound_lower_cur = Inf*ones(1, nfuns);		
				bound_upper_cur = zeros(1, nfuns);
			end
			% Compute last d+1 Lanczos approximations:
			for j = 1:d+1
				l = its_lanc+1 - (d+1) + j;
				trestk(j, :) = matfun_proj_qf(A, Omega(:, k), fmu, V(:, 1:l), T(1:l, 1:l));
			end
			% Compute upper and lower bounds:
			for j = 1:d
				boundk = c * abs(trestk(j, :) - trestk(j+1, :));
				% Upper and lower bounds for Omega(:,k)' * P_mu * Omega(:,k):
				boundk_upper = cummin(trestk(j, :) + boundk, 'reverse');
				boundk_lower = cummax(trestk(j, :) - boundk);
				% Update bounds from last few iterations:
				if bound_type == "diff"
					bound_upper_cur = min(bound_upper_cur, boundk_upper);
					bound_lower_cur = max(bound_lower_cur, boundk_lower);
				elseif bound_type == "diffsafe"
					bound_upper_cur = max(bound_upper_cur, boundk_upper);
					bound_lower_cur = min(bound_lower_cur, boundk_lower);
				end
			end
			trestk = trestk(d, :);		% for consistency
		end

		% Update averages:
		trest = trest + trestk/N_hutch;
		bound_upper = bound_upper + bound_upper_cur/N_hutch;
		bound_lower = bound_lower + bound_lower_cur/N_hutch;
	end

end