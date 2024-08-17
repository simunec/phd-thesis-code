% Test performance of gapfinder with decreasing gap size

clear;
n = 30000;

thetas = [0.1, 0.05, 0.025, 0.01, 0.005, 0.0025];
% thetas = [0.025];

thetas_actual = zeros(size(thetas));
t = zeros(size(thetas));
t_eig = zeros(size(thetas));

rng(0);
for idx_theta = 1:length(thetas)

	theta = thetas(idx_theta)
	a1 = 1;
	b1 = 1e3;


	b2 = 1e4;
	x = 2*(b2-b1)*theta/(1+theta);
	a2 = b1 + x;

	ev = [logspace(log10(a1), log10(b1), 2*floor(n/3)), logspace(log10(a2), log10(b2), n - 2*floor(n/3))].';

	v1 = randn(n, 1);
	v2 = randn(n-1, 1);

	T = spdiags(v1, 0, n, n);
	T(2:n+1:end) = v2;
	T(n+1:n+1:end) = v2;

	A = spdiags(ev, 0, n, n);
	A = A + T;

	tic;
	ee = eig(A);
	t_eig(idx_theta) = toc;

	% Check actual theta:
	l1 = min(ee);
	l2 = max(ee(ee < b1+x/2));
	l3 = min(ee(ee >= b1+x/2));
	l4 = max(ee);
	mu_actual = (l2 + l3)/2;
	theta_actual =  min(abs(mu_actual-l2), abs(mu_actual-l3)) / max(abs(mu_actual-l1), abs(mu_actual-l4));

	theta;
	thetas_actual(idx_theta) = theta_actual;

	mus = logspace(log10(a1), log10(b2), 10000); 
	delta = 1e-2;
	bound_type = "diffsafe";
	c = 2;
	d = 3;

	tic;
	[gaps, trest, trest_upper, trest_lower, ~, its_lanc(idx_theta)] = gapfinder_main(A, mus, delta, theta, bound_type, d, c);
	t(idx_theta) = toc;

	% Plots
	% figure;
	% plot(mus, trest, '-.r'); hold on;
	% plot(mus, trest_upper, '-.b'); hold on;
	% plot(mus, trest_lower, '-.b'); hold on;

	gap_true{idx_theta} = [l2, l3];
	eigcount_true(idx_theta) = sum(ee <= l2);

	% Find largest gap (required in case other smaller gaps are found):
	gapwidth = [];
	for j = 1:length(gaps)
		gapwidth(j) = gaps{j}(2) - gaps{j}(1);
	end
	[~, j] = max(gapwidth);
	idx = find(mus >= gaps{j}(1), 1);
	eigcount_est(idx_theta) = round(trest(idx));
	gaps_est{idx_theta} = gaps{j};

end


% Construct table:
% theta	| true gap | estimated gap | estimated eigcount | lanczos its | time | time_eig
table = "";
for j = 1:length(thetas)
	tablerow = sprintf("%.4f & [%.2f, %.2f] & [%.2f, %.2f] & %d & %d & %.3f & %.3f \\\\", thetas(j), gap_true{j}, gaps_est{j}, eigcount_est(j), its_lanc(j), t(j), t_eig(j));

	table = sprintf("%s\n%s", table, tablerow);
end

fname = "tables/table_gapfinder_gapscale.txt";
fid = fopen(fname, 'a');
fprintf(fid, "%s\n\n", table);
fclose(fid);




