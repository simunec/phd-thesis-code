% Generates the legend for all the figures

figure;
p1 = semilogy(NaN);	hold on;
p12 = semilogy(NaN, '-.');	hold on;
p13 = semilogy(NaN, 'x'); 	hold on;
p14 = semilogy(NaN, '+'); 	hold on;

p2 = semilogy(NaN);	hold on;
p22 = semilogy(NaN, '-.');	hold on;
p23 = semilogy(NaN, 'o');	hold on;
p24 = semilogy(NaN, '+');	hold on;

p3 = semilogy(NaN);	hold on;
p32 = semilogy(NaN, '-.');	hold on;
p33 = semilogy(NaN, 's');	hold on;
p34 = semilogy(NaN, '+');	hold on;

p4 = semilogy(NaN); 	hold on;
p42 = semilogy(NaN, '-.'); 	hold on;
p43 = semilogy(NaN, '^');	hold on;
p44 = semilogy(NaN, '+');	hold on;

pEDS = semilogy(NaN, '-.k'); 


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
axis off;

legend_handle = legend('polynomial', 'polynomial w/ shift', ...
'polynomial w/ proj.', 'polynomial w/ impl. proj.',  ...
'S\&I, $\xi = -t^{-2/\alpha}$', ...
'S\&I w/ shift, $\xi = -t^{-2/\alpha}$', ... 
'S\&I w/ proj., $\xi = -t^{-2/\alpha}$ ', ...
'S\&I w/ impl. proj., $\xi = -t^{-2/\alpha}$',  ...
'S\&I, $\xi = -\sqrt{|\lambda_{2} \lambda_{n}|}$', ...
'S\&I w/ shift, $\xi = -\sqrt{|\lambda_{2} \lambda_{n}|}$', ...
'S\&I w/ proj., $\xi = -\sqrt{|\lambda_{2} \lambda_{n}|}$', ...
'S\&I w/ impl. proj., $\xi = -\sqrt{|\lambda_{2} \lambda_{n}|}$', ...
'EDS', 'EDS w/ shift', ...
'EDS w/ proj.', 'EDS w/ impl. proj.', ...
'EDS Cauchy-Stieltjes rate', ...
'Interpreter', 'LaTeX','Location', 'eastoutside');	% 'eastoutside'

set(gcf,'Position',[0,0,1024,1024]);

set(gcf,'Position',(get(legend_handle,'Position')...
    .*[0, 0, 1, 1].*get(gcf,'Position')));
set(legend_handle,'Position',[0.01,0.01, 0.98,0.98]);
set(gcf, 'Position', get(gcf,'Position') + [500, 400, 0, 0]);


figure;
p1 = semilogy(NaN, '-.');	hold on;
%~ p12 = semilogy(nerr1s(1:len1), '-.');	hold on;
p1c = semilogy(NaN);	hold on;

p2 = semilogy(NaN, '-.');	hold on;
%~ p22 = semilogy(nerr2s(1:len2), '-.');	hold on;
p2c = semilogy(NaN);	hold on;

p3 = semilogy(NaN, '-.');	hold on;
%~ p32 = semilogy(nerr3s(1:len3), '-.');	hold on;
p3c = semilogy(NaN);	hold on;

p4 = semilogy(NaN, '-.'); 	hold on;
%~ p42 = semilogy(nerrgs(1:leng), '-.'); 	hold on;
p4c = semilogy(NaN);	hold on;


p2.Color = p1c.Color;
p4.Color = p3.Color;
p3.Color = p2c.Color;
p1c.Color = p1.Color;
p2c.Color = p2.Color;
p3c.Color = p3.Color;
p4c.Color = p4.Color;
axis off;

legend_handle = legend('polynomial, no correction', ...
'polynomial, with correction  (6.3.3)', ...
'S\&I, $\xi = -t^{-2/\alpha}$, no correction', ...
'S\&I, $\xi = -t^{-2/\alpha}$,  with correction (6.3.3)', ... 
'S\&I, $\xi = -\sqrt{|\lambda_{2} \lambda_{n}|}$, no correction', ...
'S\&I, $\xi = -\sqrt{|\lambda_{2} \lambda_{n}|}$, with correction (6.3.3)', ...
'EDS, no correction', 'EDS, with correction (6.3.3)', ...
'Interpreter', 'LaTeX','Location', 'eastoutside');	% 'eastoutside'

set(gcf,'Position',[0,0,1024,1024]);

set(gcf,'Position',(get(legend_handle,'Position')...
    .*[0, 0, 1, 1].*get(gcf,'Position')));
set(legend_handle,'Position',[0.01,0.01, 0.98,0.98]);
set(gcf, 'Position', get(gcf,'Position') + [500, 400, 0, 0]);
