% Test performance of gapfinder with increasing matrix size, with fixed number of Lanczos iterations

clear;

load linverse;

A = Problem.A;
n = size(A, 1);

rng(0);
fprintf("n = %d\n", n);

tic;
ee = eig(full(A), 'vector');
t_eig = toc;
fprintf("Eigenvalue computation time: %.2f\n", t_eig)

mus = linspace(min(ee), max(ee), 1000); 
delta = 1e-2;
its_lanc = 100;
bound_type = "diffsafe";
c = 2;
d = 3;

rng(0);
tic;
[gaps, trest, trest_upper, trest_lower, ~] = gapfinder_main(A, mus, delta, its_lanc, bound_type, d, c);
t_diff = toc;

bound_type = "residue";
rng(0);
tic;
[gaps_res, trest_res, trest_upper_res, trest_lower_res, ~] = gapfinder_main(A, mus, delta, its_lanc, bound_type, d, c);
t_res = toc;


% Find truegaps corresponding to found gaps:
for j = 1:length(gaps)
	mu = 0.5*(gaps{j}(2) + gaps{j}(1));
	idx = find(ee > mu, 1);
	truegaps{j} = [ee(idx-1), ee(idx)];

end

% Plots:
figure;
p1 = plot(ee, 'o'); 
p1.MarkerSize = 3;
xlim([-400, length(ee) + 400]);
ylim([min(ee)-1, max(ee)+1]);
yticks([-5, 0, 5, 10, 15]);

set(gcf,'PaperPositionMode','auto');
set(gcf,'PaperSize', [6 4]);
set(gcf, "PaperPosition", [0 0 6 4]);

fname = "plots/realmat_eig";
print(fname, '-depsc2');


figure;
plot(mus, trest, '-.r'); hold on;
plot(mus, trest_upper, '-.b'); hold on;
plot(mus, trest_lower, '-.b'); hold on;
xlim([min(mus) - 1, max(mus) + 1]);
ylim([min(trest) - 600, max(trest) + 600]);
xlabel("\mu");

set(gcf,'PaperPositionMode','auto');
set(gcf,'PaperSize', [6 4]);
set(gcf, "PaperPosition", [0 0 6 4]);

fname = "plots/gapfinder_realmat_diff";
print(fname, '-depsc2');


figure;
plot(mus, trest_res, '-.r'); hold on;
plot(mus, trest_upper_res, '-.b'); hold on;
plot(mus, trest_lower_res, '-.b'); hold on;
xlim([min(mus) - 1, max(mus) + 1]);
ylim([min(trest) - 600, max(trest) + 600]);
xlabel("\mu");

set(gcf,'PaperPositionMode','auto');
set(gcf,'PaperSize', [6 4]);
set(gcf, "PaperPosition", [0 0 6 4]);

fname = "plots/gapfinder_realmat_res";
print(fname, '-depsc2');

% Construct table:
% first gap | second gap | third gap | time
tablerow1 = sprintf("eig & [%.3f, %.3f] & [%.3f, %.3f] & [%.3f, %.3f] & %.3f \\\\", truegaps{1:3}, t_eig);
tablerow2 = sprintf("consec.~diff. & [%.3f, %.3f] & [%.3f, %.3f] & [%.3f, %.3f] & %.3f \\\\", gaps{1:3}, t_diff);
tablerow3 = sprintf("residues & [%.3f, %.3f] & [%.3f, %.3f] & [%.3f, %.3f] & %.3f \\\\", gaps_res{1:3}, t_res);
table = sprintf("%s\n%s\n%s", tablerow1, tablerow2, tablerow3);

fname = "tables/table_gapfinder_realmat.txt";
fid = fopen(fname, 'a');
fprintf(fid, "%s\n\n", table);
fclose(fid);
