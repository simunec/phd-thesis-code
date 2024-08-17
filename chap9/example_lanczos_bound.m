% Test Lanczos approximation of a quadratic form and compare with a priori and a posteriori error bounds

clear;
n = 1000;

a1 = 0;
b1 = a1 + 10;

a2 = b1 + 1;
b2 = a2 + 9;

a3 = b2 + 5;
b3 = a3 + 15;

n1 = floor(n/3);
n2 = n1;
n3 = n - 2*n1;

ev = [chebpts(n1, [a1, b1]); chebpts(n2, [a2, b2]); chebpts(n3, [a3, b3])];
% ev = [logspace(log10(a1+1), log10(b1), n1),  logspace(log10(a2), log10(b2), n2), logspace(log10(a3), log10(b3), n3)];
lmin = min(ev);
lmax = max(ev);

profile on;

rng(0);
[A, U, D] = generate_testmatrix(ev);		% A = U*D*U'
b = randn(n, 1);
b = b / norm(b);

% mu = 22;
% its_lanc = 100;

% % Second test:
mu = 10.6;
its_lanc = 400;	


offset = 1;		% offset for consecutive difference error estimate
fmu = @(z) real(z) <= mu;
dfmu = @(z) 0;


qf_exact = (b'*U) * diag(fmu(diag(D))) * (U'*b);

[V, T] = lanczos(A, b, its_lanc);
qf = zeros(1, its_lanc);
qf_bound_post = zeros(1, its_lanc);
qf_est_diff = zeros(1, its_lanc);
for j = 1:its_lanc	
	% Approximation with Lanczos:
	[qf(j), UU, DD] = matfun_proj_qf(A, b, fmu, V(:, 1:j), T(1:j, 1:j));	
	
	% A posteriori error bound:
	[~, ~, ~, theta, a2b2, abc2] = cauchy_integral_residues_setup(UU, DD, T(j+1, 1:j));
	lambda = linspace(min(theta), max(theta), 10000);
	fmus = {fmu}; 	dfmus = {dfmu};		% for compatibility with next function
	[qf_bound_post(j), ~] = cauchy_integral_residues_quadform(a2b2, abc2, theta, fmus, dfmus, lambda);

	% Consecutive difference estimate:
	if j > offset
		qf_est_diff(j-offset) = abs(qf(j-offset) - qf(j));
	end
end
qf_est_diff(its_lanc-offset+1:its_lanc) = inf;
qf_bound_post = qf_bound_post .* norm(b)^2;		% multiply by norm

% A priori error bound:
mm = 2:its_lanc;
% Compute theta, relative gap width
l1 = min(ev);
l2 = max(ev(ev < mu));
l3 = min(ev(ev >= mu));
l4 = max(ev);
mu_avg = (l2 + l3)/2;
theta =  min(abs(mu-l2), abs(mu-l3)) / max(abs(mu-l1), abs(mu-l4));
C = (1-theta)/sqrt(pi*theta) + 1;
theta_avg =  min(abs(mu_avg-l2), abs(mu_avg-l3)) / max(abs(mu_avg-l1), abs(mu_avg-l4));
C_avg = (1-theta_avg)/sqrt(pi*theta_avg) + 1;
qf_bound_prior(1) = inf;
qf_bound_prior(2:its_lanc) = (C*norm(b)^2 ./ sqrt(mm - 1)) .* ((1-theta)/(1+theta)).^(mm-1);
qf_bound_prior_avg(1) = inf;
qf_bound_prior_avg(2:its_lanc) = (C_avg*norm(b)^2 ./ sqrt(mm - 1)) .* ((1-theta_avg)/(1+theta_avg)).^(mm-1);

profile viewer;

%* Plots:
figure;
leg = strings(0);
semilogy(abs(qf_exact - qf), '-k');	hold on;
leg(end+1) = "error";
semilogy(qf_bound_post, '-.b');		hold on;
leg(end+1) = "a posteriori bound";
semilogy(qf_est_diff, '-.g');	hold on;
leg(end+1) = "consec.~diff.~estimate";
semilogy(qf_bound_prior, '--r');	hold on;
leg(end+1) = "a priori bound";
semilogy(qf_bound_prior_avg, '--m');	hold on;
leg(end+1) = "a priori bound, central $\mu$";

legend(leg, 'interpreter', 'latex', 'location', 'northeast');

ylim([1e-15, 1e5]);
xlim([-5, its_lanc+5]);
xlabel("$m$", 'interpreter', 'latex');
ylabel("error", 'interpreter', 'latex');
set(gcf,'PaperPositionMode','auto');
set(gcf,'PaperSize', [6 4]);
set(gcf, "PaperPosition", [0 0 6 4]);

fname = sprintf("plots/example_lanczos_bound%d", its_lanc);
print(fname, '-depsc2');
