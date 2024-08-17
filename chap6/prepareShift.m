function [x,z] = prepareShift(L, type)
	% [x,z] = prepareShift(L, type) 
	% Computes vectors x and z such that L + z*x' is nonsingular, L z = 0 and x' L = 0
	% The vectors x and z are nonnegative and are normalized so that x' z = 1.
	% L should be the transpose of the graph Laplacian of a strongly connected graph. 
	% If type = 'u', the graph is undirected and x = z = ones(n) up to normalization.
	% Otherwise, x = ones(n) and z is computed with an LDU factorization of L.
	% z is always normalized so that sum(z) = 1.
	if nargin < 2
		type = "";
	end
	
	n = size(L, 1);
	x = ones(n,1);
	if type == "u"
		z = ones(n,1)./n;
	else	
		perm = amd(L);
		[l,u] = lu(L(perm,perm));
		uu = u(1:end-1, 1:end-1);
		z = uu \ u(1:end-1, end);
		z = [-z; 1];
		z(perm) = z./sum(z);
	end
end
