function [A, graphname] = extractLCC(graphname, type, pathname)
	% A = extractLCC(graphname, type, pathname)
	% Returns the adjacency matrix A of the largest connected component of the graph "graphname.mat".
	% optional argument pathname can be used ignore other arguments
	% and instead do load(pathname)

	% convert graphname to char vector
	graphname = char(graphname);
	% Find the file with the graph
	if nargin < 2
		type = "";
	end
	if nargin < 3
		% assume that graph is already in the path
		pathname = graphname;
	end
	load(pathname);
	
	% Load the adjacency matrix of graphname
	A = Problem.A;	
	
	% Extract the largest connected component
	if type == "u"
		G = graph(A);			
	else 
		G = digraph(A);
	end
	[bins, binsize] = conncomp(G);	% strongly connected components
	idx = find(bins == mode(bins));	% indices of largest connected component
	GG = subgraph(G, idx);
	GG = simplify(GG);		% remove multiple edges
	A = adjacency(GG);
end
