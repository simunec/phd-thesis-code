function [gaps, idx_widest] = detect_gaps(mu, bound_upper, bound_lower, tol)
	% [gaps, idx_widest] = detect_gaps(mu, bound_upper, bound_lower, tol)
	% Finds gaps in the spectrum of a matrix, using upper and lower bounds
	% on the approximation of trace(P_mu) using a single Hutchinson Gaussian vector
	% Input:
	%	mu				vector of shift parameters for spectral projector P_mu
	%	bound_upper		upper bound on trace^H(P_mu)
	%	bound_lower		lower bound on trace^H(P_mu)
	%	tol				tolerance for detecting gaps
	% Output:
	%	gaps			cell array with pairs [a_i, b_i] (extrema of found gaps)
	%	idx_widest		vector with i-th entry equal to the index of the largest possible gap starting from mu(i)
	%
	% The gaps are found according to the following criteria:
	%	1. the union of all intervals [mu(i), mu(i+1)] such that bound_upper(i+1) - bound_lower(i) < tol
	%	has a probability smaller than delta to contain an eigenvalue of A 
	%	(this corresponds to a small jump in trace^H(P_mu)) 
	%	2. however this is not enough to determine a gap in the spectrum with probability 1-delta, because an
	%	interval without eigenvalues must contain a horizontal line of trace^H(P_mu)
	%	3. so for each i we also need to find the largest index idx_widest(i) such that the interval
	%	[mu(i), mu(idx_widest(i))] can contain a horizontal line, i.e. such that
	%	bound_upper(i) > bound_lower(idx_widest(i))
	%	4. if an interval satisfies (1) but not (3), then it must certainly contain an eigenvalue of A,
	%	although this situation is very unlikely to happen in practice

	N = length(mu);
	mu = mu(:)';
	bound_upper = bound_upper(:)';
	bound_lower = bound_lower(:)';
	% correct numerical errors so that bound_upper(i) >= bound_lower(i) for all i:
	bound_upper = max(bound_upper, bound_lower);
	idx_widest = zeros(size(mu));
	gaps = cell(0);

	% Find largest possible horizontal line of trace^H(P_mu) starting from each point in mu:
	idx_u = 1;	idx_l = 1;
	% this costs O(N), since in each iteration we increase either idx_u or idx_l
	while idx_l <= N && idx_u <= N
		if bound_upper(idx_u) >= bound_lower(idx_l) - 1e-10		
			% above comparison may need some looseness tolerance for oscillations due to machine precision
			% we can draw a horizontal line from idx_u to idx_l
			% (this always happens if idx_u == idx_l)
			% idx_widest(idx_u) = idx_l;
			if idx_l == N
				% Reached final index, update all remaining indices up to N
				idx_widest(idx_u:N) = idx_l;
			end
			idx_l = idx_l+1;
		else
			% we cannot draw a horizontal line from idx_u to idx_l
			% go to idx_u+1, noticing that we can still do a horizontal line from idx_u+1 to idx_l-1
			idx_widest(idx_u) = idx_l-1;
			idx_u = idx_u+1;
		end
	end

	% Find consecutive intervals without eigenvalues:
	gapmask = abs(bound_upper(2:N) - bound_lower(1:N-1)) < tol;
	% true at indices where there should be no eigenvalues, false otherwise
	idx_mask = 1;
	foundfirst = false;		% true if we have found a streak of true inside gapmask
	while idx_mask <= N-1
		if foundfirst == false
			% look for a true value in gapmask and use it to set max index that allows horizontal lines
			if gapmask(idx_mask) == false
				idx_mask = idx_mask+1;
			else
				% found starting index for a gap!
				foundfirst = true;
				idx_first = idx_mask;				% gapmask(i) corresponds to interval [mu(i), mu(i+1)]
				max_idx = idx_widest(idx_mask)-1;	% horizontal line until mu(j) means until the interval gapmask(j-1)
				% it may happen here that max_idx = idx_mask - 1, when it is impossible to draw a horizontal line
				% in this case, the interval [mu(idx_mask), mu(idx_mask+1)] cannot be a gap, and it probably
				% contains an eigenvalue with a small jump associated to it (unlikely scenario)
				if max_idx < idx_mask
					% discard potential gap and go back to looking for starting index
					foundfirst = false;
				end
				idx_mask = idx_mask+1;
			end
		else
			% keep increasing idx_mask until either gapmask is false or idx_mask exceeds max_idx
			if gapmask(idx_mask) == true && idx_mask <= max_idx
				% add current gap if reached end of gapmask
				if idx_mask == N-1
					gaps{end+1} = [mu(idx_first), mu(idx_mask+1)];
				end
				idx_mask = idx_mask+1;
			else
				% cannot further increase current gap
				gaps{end+1} = [mu(idx_first), mu(idx_mask)];

				% note that here idx_mask comes from (idx_mask - 1) + 1
				% (previous index increased by 1 because gapmask(i) corresponds to interval [mu(i), mu(i+1)])
				foundfirst = false;
				if gapmask(idx_mask) == false
					% should not increase if true, since the "next" gap can start right away
					idx_mask = idx_mask+1;
				end
			end
		end
	end
		

end