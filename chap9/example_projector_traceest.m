% Test basic Hutchinson + Lanczos trace estimator on a small problem

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

mus = linspace(lmin-5, lmax+5, 2000); 
% N_hutch = 1;		% Hutchinson vectors
N_hutch = 10;		% Hutchinson vectors
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
krylovfun = @(A, b) lanczos(A, b, its_lanc);
trest_hl_avg = zeros(1, nfuns);

for k = 1:N_hutch
	fprintf("Hutchinson + Lanczos %d\n", k);

	%* Compute Lanczos approximation:
	[V, T] = krylovfun(A, Omega(:, k));
	[trest_hl, UU, DD] = matfun_proj_qf(A, Omega(:, k), fmu, V(:, 1:its_lanc), T(1:its_lanc, 1:its_lanc));	

	%* Update average:
	trest_hl_avg = trest_hl_avg + trest_hl/N_hutch;
end

%* Plots:
figure;
leg = strings(0);
plot(mus, n_mu, '-k');	hold on;
leg(end+1) = "$\mathrm{n}_\mathrm{e}(\mu)$";
plot(mus, trest_h, 'r');	hold on;
leg(end+1) = sprintf("Hutchinson $s=%d$, exact $\\textit{\\textbf{x}}^T P_{\\mu} \\textit{\\textbf{x}}$", N_hutch);
phl_avg = plot(mus, trest_hl_avg, 'b');	hold on;
leg(end+1) = sprintf("Hutchinson $s=%d$, Lanczos $m=%d$", N_hutch, its_lanc);
legend(leg, 'interpreter', 'latex', 'location', 'southeast');

ylim([-20, 670]);
xlim([-5, 65]);
xlabel("\mu");
set(gcf,'PaperPositionMode','auto');
set(gcf,'PaperSize', [6 4]);
set(gcf, "PaperPosition", [0 0 6 4]);


fname = sprintf("plots/example_projector_traceest%d", N_hutch);
print(fname, '-depsc2');
