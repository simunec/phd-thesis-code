% Test effect of correction on standard Krylov methods

% PARAMETERS: 
	graphname = "Roget";
	type = "d";
	t = 1;
	alpha = 0.5;
	truesol_method = 'eig';			% uses eig function
	maxit1 = 150;					% polynomial
	maxit2 = 60;					% S&I Moret-Novati
	maxit3 = 60;					% S&I spectral
	maxitg = 40;					% rational EDS

	fv = @(B, w) fracexpmvec(B, alpha, t, w);		% exp(-t B^alpha)*w

	tic;
	A = extractLCC(graphname, type);
	n = size(A,1);		
	L = spdiags(A*ones(n,1), 0, n, n) - A;
	L = sparse(L');	
	temp = toc;
	fprintf("Graph preprocessing took %.2f seconds\n\n", temp);	
	% NOTE: L is the transpose of the graph Laplacian!

	% random initial probability vector
	v = rand(n,1);
	v = v./sum(v);	

	%%% ------ FOR COMPUTING SHIFT ONLY ONCE ACROSS ALL METHODS ------ %%%
	% needed by correction function
	[x, z] = prepareShift(L, type);
	% x = ones(n,1) and sum(z) = 1, L z = 0 (L is transposed Laplacian)

	% L_rk1shifted = L + theta*z*x';	% exact shifted Laplacian 
	% f(L)v = f(L + theta zx')v - [f(theta)-f(0)]z x'v
	% backshift is such that f(L) v = f(L_shifted) v - backshift
	% backshift = (f(theta)-f(0))*z*(x'*v);
	%%% -------------------------------------------------------------- %%%

	% Compute spectral choice for shift delta, used later:
	tic;
	lambda_n = eigs(L, 1);	% largest modulus eigenvalue
	lambda_2 = eigs(L + speye(n), 2, 'smallestabs');
	% second smallest eigenvalue of L + I(the other one should be zero)
	lambda_2 = max(lambda_2) - 1;	% now of L
	temp = toc;
	fprintf("Computing extremal eigenvalues took %.2f seconds\n\n", temp);

	delta_spec = sqrt(abs(lambda_2)*abs(lambda_n));

	switch truesol_method
		case 'eig'	
		%%% ---------- Compute true solution for graph Laplacians:
			fprintf("Computing actual true solution using eig: ");
			tic;
			[Q,D] = eig(full(L));
			if type == "u"
				y_true =  Q*diag( exp(-t*diag(D).^(alpha)) ) * (Q'*v);
			else
				y_true = Q*diag( exp(-t*diag(D).^(alpha)) ) * (Q\v);
			end
			y_true = real(y_true);
			y_true = correction(y_true, z);
			t_true = toc;
			fprintf("%.2f seconds\n", t_true);
		case 'projection'
			fprintf("Computing 'true' solution using EDS RK method with implicit projection: ");
			its = 60;
			tic;
			poles_true = eds_cauchy_st(0.99*lambda_2, 1.01*lambda_n, its, 1/sqrt(3));	
			%%% ---------- IMPLICIT PROJECTION CHOICE: 
			w = v - z;
			fz = fv(0,1)*z; 
			y_ref = runRatKrylov(L, w, poles_true, fv, "");
			checknorms = zeros(size(y_ref,2)-1, 1);
			for k = 1:length(checknorms)
				checknorms(k) = norm(y_ref(:,k+1) - y_ref(:,k));
			end
			[~, idx] = min(checknorms);
			y_true = y_ref(:, idx) + fz;		% hopefully close to true solution
			t_true = toc;
			fprintf("%.2f seconds\n", t_true);
			% figure;
			% semilogy(checknorms, 'k');
		otherwise
			disp("Unknown truesol_method");
			exit();
	end

	%%% ----- Standard methods ----- %%%
	disp("Starting standard methods");
	% Polynomial: all poles at infinity
	poles = Inf*ones(1, maxit1);
	[y1, t1] = runRatKrylov(L, v, poles, fv, "", maxit1);
	nerr1 = getErrorNorms(y_true, y1);
	% With correction:
	y1c = correction(y1, z);
	nerr1c = getErrorNorms(y_true, y1c);

	% Shift-and-Invert: single repeated pole
	delta = t^(-2/alpha); % MORET-NOVATI CHOICE
	poles = -delta*ones(1, maxit2);
	[y2, t2] = runRatKrylov(L, v, poles, fv, "", maxit2);
	nerr2 = getErrorNorms(y_true, y2);
	% With correction:
	y2c = correction(y2, z);
	nerr2c = getErrorNorms(y_true, y2c);

	delta = delta_spec; % SPECTRAL CHOICE
	poles = -delta*ones(1, maxit3);
	[y3, t3] = runRatKrylov(L, v, poles, fv, "", maxit3);
	nerr3 = getErrorNorms(y_true, y3);
	% With correction:
	y3c = correction(y3, z);
	nerr3c = getErrorNorms(y_true, y3c);

	% General Krylov with asymptotically optimal poles (EDS)
	poles_eds = eds_cauchy_st(0.99*lambda_2, 1.01*lambda_n, maxitg+1);
	poles = poles_eds(1:maxitg);
	[yg, tg] = runRatKrylov(L, v, poles, fv, "", maxitg);
	nerrg = getErrorNorms(y_true, yg);
	% With correction:
	ygc = correction(yg, z);
	nerrgc = getErrorNorms(y_true, ygc);
	
	% Plot results:
	len1 = length(nerr1);		% polynomial
	len2 = length(nerr2);		% S&I Moret-Novati
	len3 = length(nerr3);		% S&I spectral
	leng = length(nerrg);		% rational Cauchy-Stieltjes

	figure;
	p1 = semilogy(nerr1(1:len1), '-.');		hold on;
	p1c = semilogy(nerr1c(1:len1));			hold on;

	p2 = semilogy(nerr2(1:len2), '-.');		hold on;
	p2c = semilogy(nerr2c(1:len2));			hold on;

	p3 = semilogy(nerr3(1:len3), '-.');		hold on;
	p3c = semilogy(nerr3c(1:len3));			hold on;

	p4 = semilogy(nerrg(1:leng), '-.'); 	hold on;
	p4c = semilogy(nerrgc(1:leng)); 		hold on;

	p4.Color = p3.Color;
	p4c.Color = p3.Color;
	p3.Color = p2c.Color;
	p3c.Color = p2c.Color;
	p2.Color = p1c.Color;
	p2c.Color = p1c.Color;
	p1.Color = p1.Color;
	p1c.Color = p1.Color;

	xlabel('iterations', 'interpreter', 'latex');
	ylabel('relative error', 'interpreter', 'latex');
	ylim([0.5e-14, 2e0]);
	xlim([-5, max([len1, len2, len3, leng]) + 5]);

	title(sprintf('%s -- $n = %d$, $t = %.1f$, $\\alpha = %.2f$', graphname, n, t, alpha),'Interpreter','LaTeX')

	set(gcf,'PaperPositionMode','auto');
	set(gcf,'PaperSize', [6 4]);
	set(gcf, "PaperPosition", [0 0 6 4]);

	figurename = sprintf("plots/%s_%.2f_%.2f.eps", graphname, alpha, t);
	print(figurename, '-depsc2');
		
	function y = correction(y_temp, z)
		% y = correction(y_temp, z)
		% removes error along the direction z, where L z = 0
		% works assuming sum(z) = 1
		% the correction is computed by enforcing sum(y(:, k)) = 1 for all k
		% has good effect experimentally if no shift or projection is used 
		corr = (sum(y_temp,1) - 1);		% row vector of corrections
		y = y_temp - z*corr;			% corrected vectors
	end
	
	function nerr = getErrorNorms(y_true, y)
		% nerr = getErrorNorms(y_true, y)
		% finds the relative error norms: 
		% norm(y_true - y(:,k))/norm(y_true) for all k
		err = y - y_true;
		nerr = zeros(size(err,2), 1);
		ny_true = norm(y_true);
		for k = 1:length(nerr)
			nerr(k) = norm(err(:,k))/ny_true;
		end
	end
	