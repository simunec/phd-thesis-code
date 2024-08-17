% Plot results


% Cauchy-Stieltjes convergence rate
EDSrate_cs = exp(-pi^2 / log(4*4*lambda_n/lambda_2));

%%% ---------- PLOTTING ---------- %%%
% assume that all curves related to the same pole choice have the same length:
len1 = length(nerr1);		% polynomial
len2 = length(nerr2);		% S&I Moret-Novati
len3 = length(nerr3);		% S&I spectral
leng = length(nerrg);		% rational Cauchy-Stieltjes

figure;
p1 = semilogy(nerr1(1:len1));	hold on;
p12 = semilogy(nerr1s(1:len1), '-.');	hold on;
p13 = semilogy([8:8:len1], nerr1p(8:8:len1), 'x'); 	hold on;
p14 = semilogy([4:8:len1], nerr1ip(4:8:len1), '+'); 	hold on;

p2 = semilogy(nerr2(1:len2));	hold on;
p22 = semilogy(nerr2s(1:len2), '-.');	hold on;
p23 = semilogy([8:8:len2], nerr2p(8:8:len2), 'o');	hold on;
p24 = semilogy([4:8:len2], nerr2ip(4:8:len2), '+');	hold on;

p3 = semilogy(nerr3(1:len3));	hold on;
p32 = semilogy(nerr3s(1:len3), '-.');	hold on;
p33 = semilogy([8:8:len3], nerr3p(8:8:len3), 's');	hold on;
p34 = semilogy([4:8:len3], nerr3ip(4:8:len3), '+');	hold on;

p4 = semilogy(nerrg(1:leng)); 	hold on;
p42 = semilogy(nerrgs(1:leng), '-.'); 	hold on;
p43 = semilogy([8:8:leng], nerrgp(8:8:leng), '^');	hold on;
p44 = semilogy([4:8:leng], nerrgip(4:8:leng), '+');	hold on;

pEDS_cs = semilogy(EDSrate_cs.^[1:30], '-.k'); 

p4.Color = p2.Color;
p3.Color = p14.Color;
p2.Color = p12.Color;

p12.Color = p1.Color;
p13.Color = p1.Color;
p14.Color = p1.Color;
p22.Color = p2.Color;
p23.Color = p2.Color;
p24.Color = p2.Color;
p32.Color = p3.Color;
p33.Color = p3.Color;
p34.Color = p3.Color;
p42.Color = p4.Color;
p43.Color = p4.Color;
p44.Color = p4.Color;


xlabel('iterations', 'interpreter', 'latex');
ylabel('relative error', 'interpreter', 'latex');
ylim([0.5e-15, 2e0]);
xlim([-5, max([len1, len2, len3, leng]) + 5]);

title(sprintf('%s -- $n = %d$, $t = %.1f$, $\\alpha = %.2f$', graphname, n, t, alpha),'Interpreter','LaTeX')

set(gcf,'PaperPositionMode','auto');
set(gcf,'PaperSize', [6 4]);
set(gcf, "PaperPosition", [0 0 6 4]);

figurename = sprintf("plots/%s_%.2f_%.2f.eps", graphname, alpha, t);
print(figurename, '-depsc2');
