% Test performance of gapfinder with increasing matrix size, with fixed number of Lanczos iterations

nn = [5000, 10000, 20000, 40000, 80000];

% nn = 5000;

theta = 0.01;

t = zeros(size(nn));
t_eig = zeros(size(nn));

rng(0);
for idx_n = 1:length(nn)

	n = nn(idx_n);
	fprintf("n = %d\n", n);

	a1 = 1;
	b1 = 1e3;

	b2 = 1e4;
	x = 2*(b2-b1)*theta/(1+theta);
	a2 = b1 + x;

	ev = [logspace(log10(a1), log10(b1), floor(n/2)), logspace(log10(a2), log10(b2), n - floor(n/2))].';

	v1 = randn(n, 1);
	v2 = randn(n-1, 1);

	T = spdiags(v1, 0, n, n);
	T(2:n+1:end) = v2;
	T(n+1:n+1:end) = v2;

	A = spdiags(ev, 0, n, n);
	A = A + T;

	tic;
	ee = eig(A);
	t_eig(idx_n) = toc;

	% Check actual theta:
	l1 = min(ee);
	l2 = max(ee(ee < b1+x/2));
	l3 = min(ee(ee >= b1+x/2));
	l4 = max(ee);
	mu_actual = (l2 + l3)/2;
	theta_actual =  min(abs(mu_actual-l2), abs(mu_actual-l3)) / max(abs(mu_actual-l1), abs(mu_actual-l4));

	mus = logspace(log10(a1), log10(b2), 10000); 
	delta = 0.5;
	its_lanc = 250;
	bound_type = "diffsafe";
	c = 2;
	d = 3;

	tic;
	[gaps, trest, trest_upper, trest_lower, gaps_rough] = gapfinder_main(A, mus, delta, its_lanc, bound_type, d, c);
	t(idx_n) = toc;

	% Plots
	% figure;
	% plot(mus, trest, '-.r'); hold on;
	% plot(mus, trest_upper, '-.b'); hold on;
	% plot(mus, trest_lower, '-.b'); hold on;

	gap_true{idx_n} = [l2, l3];
	eigcount_true(idx_n) = sum(ee <= l2);

	disp(length(gaps));
	% Find largest gap (required in case other smaller gaps are found):
	gapwidth = [];
	for j = 1:length(gaps)
		gapwidth(j) = gaps{j}(2) - gaps{j}(1);
	end
	[~, j] = max(gapwidth);
	idx = find(mus >= gaps{j}(1), 1);
	eigcount_est(idx_n) = round(trest(idx));
	gap_est{idx_n} = gaps{j};

end

% Construct table:
% n	| true gap | estimated gap | estimated eigcount | lanczos its | time | time eig
table = "";
for j = 1:length(nn)
	tablerow = sprintf("%d & [%.2f, %.2f] & [%.2f, %.2f] & %d & %d & %.3f & %.3f \\\\", nn(j), gap_true{j}, gap_est{j}, eigcount_est(j), its_lanc, t(j), t_eig(j));

	table = sprintf("%s\n%s", table, tablerow);
end

fname = "tables/table_gapfinder_size_fixiter.txt";
fid = fopen(fname, 'a');
fprintf(fid, "%s\n\n", table);
fclose(fid);




