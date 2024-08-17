% Plots Lanczos trace estimates with upper and lower bounds for many different values of mu
% Same setup as example_projector_traceest.m

n = 600;

a1 = 0;
b1 = a1 + 20;

a2 = b1 + 1;
b2 = a2 + 9;

a3 = b2 + 2;
b3 = a3 + 8;

a4 = b3 + 4;
b4 = a4 + 16;

n1 = floor(n/4);
n2 = n1;
n3 = n1;
n4 = n - 3*n1;

ev = [chebpts(n1, [a1, b1]); chebpts(n2, [a2, b2]); chebpts(n3, [a3, b3]); chebpts(n4, [a4, b4])];
lmin = min(ev);
lmax = max(ev);

rng(0);
[A, U, D] = generate_testmatrix(ev);		% A = U*D*U'
b = randn(n, 1);

mus = [linspace(lmin-5, lmax+5, 4000)]; 
N_hutch = 1;		% Hutchinson vectors
its_lanc = 50;		% Krylov iterations

%* Compute exact number of eigenvalues below each mu:
n_mu = zeros(1, length(mus));
for j = 1:length(mus)
	n_mu(j) = sum(diag(D) - mus(j) < 0);
end

nfuns = length(mus);
fmu = cell(nfuns, 1);
dfmu = cell(nfuns, 1);
for k = 1:nfuns
	fmu{k} = @(z) (real(z) <= mus(k));		% step function
	dfmu{k} = @(z) 0;
end

%* Approximation using Hutchinson and exact matvecs with P_mu:
trest_h = zeros(1, length(mus));
trest_h_last = zeros(1, length(mus));
Omega = randn(n, N_hutch);				% Gaussian
UOmega = U'*Omega;
for j = 1:nfuns
	Dmu = diag(fmu{j}(diag(D)));
	trest_h(j) = trace(UOmega' * Dmu * UOmega)/N_hutch;
	trest_h_last(j) = trace(UOmega(:, end)' * Dmu * UOmega(:, end));
end

%* Approximation using Hutchinson + Lanczos:
% Compute Lanczos approximation:
[V, T] = lanczos(A, Omega, its_lanc+1);
d = 3;
c = 2;
trest_hl = zeros(d+1, nfuns);
for j = 1:d+1
	l = its_lanc+1 - (d+1) + j;
	trest_hl(j, :) = matfun_proj_qf(A, Omega, fmu, V(:, 1:l), T(1:l, 1:l));
end

% Compute trace approximation bound (consecutive differences):
bound_upper_best = Inf*ones(1, nfuns);
bound_lower_best = zeros(1, nfuns);
bound_upper_safe = zeros(1, nfuns);
bound_lower_safe = Inf*ones(1, nfuns);
for j = 1:d
	bound_err = c * abs(trest_hl(j, :) - trest_hl(j+1, :));
	% Upper and lower bounds for Omega(:,k)' * P_mu * Omega(:,k):
	bound_upper = trest_hl(j, :) + bound_err;
	bound_lower = trest_hl(j, :) - bound_err;
	bound_upper_monotone = cummin(bound_upper, 'reverse');
	bound_lower_monotone = cummax(bound_lower);
	bound_upper_best = min(bound_upper_best, bound_upper_monotone);
	bound_lower_best = max(bound_lower_best, bound_lower_monotone);

	bound_upper_safe = max(bound_upper_safe, bound_upper_monotone);
	bound_lower_safe = min(bound_lower_safe, bound_lower_monotone);
end

% Compute trace approximation bound (a posteriori bound):
bound_upper_best2 = Inf*ones(1, nfuns);
bound_lower_best2 = zeros(1, nfuns);
for il = its_lanc-d+1:its_lanc
	[trest_cur, UU, DD] = matfun_proj_qf(A, Omega, fmu, V(:, 1:il), T(1:il, 1:il));	
	[~, ~, ~, theta, a2b2, abc2] = cauchy_integral_residues_setup(UU, DD, T(il+1, 1:il));
	lambda = linspace(min(theta), max(theta), 1000);
	[bound_err2, ~] = cauchy_integral_residues_quadform(a2b2, abc2, theta, fmu, dfmu, lambda);
	bound_err2 = bound_err2 * norm(Omega)^2; % multiply by norm(Omega)^2 ~ n

	% Compute upper and lower bounds for Omega' * P_mu * Omega:
	bound_upper2 = trest_cur + bound_err2;
	bound_lower2 = trest_cur - bound_err2;
	bound_upper_monotone2 = cummin(bound_upper2, 'reverse');
	bound_lower_monotone2 = cummax(bound_lower2);

	% Update bounds from last few iterations:
	bound_upper_best2 = min(bound_upper_best2, bound_upper_monotone2);
	bound_lower_best2 = max(bound_lower_best2, bound_lower_monotone2);
end


%* Plots:
figure;
leg = strings(0);
% plot(mus, n_mu, '-k');	hold on;
% leg(end+1) = "$\mathrm{n}_\mathrm{e}(\mu)$";
plot(mus, trest_h, 'r');	hold on;
leg(end+1) = sprintf("exact $\\textit{\\textbf{x}}^T P_{\\mu} \\textit{\\textbf{x}}$");

% plot(mus, bound_upper, '--m'); hold on;
% plot(mus, bound_lower, '--m'); hold on;
plot(mus, bound_upper_best, '-.m'); hold on;
plot(mus, bound_lower_best, '-.m'); hold on;

leg(end+1) = "$\widehat{\texttt{q}}_{\mathcal{M}}^{\mathrm{U}}(\mu)$ and $\widehat{\texttt{q}}_{\mathcal{M}}^{\mathrm{L}}(\mu)$ w/ consec. diff.";
leg(end+1) = "";

plot(mus, bound_upper_safe, '-.b'); hold on;
plot(mus, bound_lower_safe, '-.b'); hold on;

leg(end+1) = "$\widehat{\texttt{q}}_{\mathcal{M}}^{\mathrm{U} \star}(\mu)$ and $\widehat{\texttt{q}}_{\mathcal{M}}^{\mathrm{L} \star}(\mu)$ w/ consec. diff.";
leg(end+1) = "";

% plot(mus, bound_upper2, '--g'); hold on;
% plot(mus, bound_lower2, '--g'); hold on;
plot(mus, bound_upper_best2, '-.g'); hold on;
plot(mus, bound_lower_best2, '-.g'); hold on;

leg(end+1) = "$\widehat{\texttt{q}}_{\mathcal{M}}^{\mathrm{U}}(\mu)$ and $\widehat{\texttt{q}}_{\mathcal{M}}^{\mathrm{L}}(\mu)$ w/ Proposition~9.2.7";
leg(end+1) = "";

legend(leg, 'location', 'southeast', 'interpreter', 'latex');

ylim([-20, 620]);
xlim([-5, 65]);

xlabel("\mu");
set(gcf,'PaperPositionMode','auto');
set(gcf,'PaperSize', [6 4]);
set(gcf, "PaperPosition", [0 0 6 4]);

fname = sprintf("plots/example_trace_bound");
print(fname, '-depsc2');

xlim([38.7, 45.3]);
ylim([388, 437]);
xticks([39:45]);
yticks([390:5:435]);
set(gcf,'PaperPositionMode','auto');
set(gcf,'PaperSize', [6 4]);
set(gcf, "PaperPosition", [0 0 6 4]);
fname2 = sprintf("plots/example_trace_bound_zoom");
print(fname2, '-depsc2');

sum(bound_upper_best - trest_h < 0)
sum(bound_upper_safe - trest_h < 0)
sum(bound_lower_best - trest_h > 0)
sum(bound_lower_safe - trest_h > 0)


