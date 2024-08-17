% testing script for KM convergence for y = exp(-t L^alpha) v
% for directed graphs
function testingScript(testcase)
% runs test corresponding to test number testcase

	% testcase 1-3: minnesota, small undirected graph
	% testcase 4-6: wiki-Vote, small directed graph
	% testcase 7-9: large undirected graphs
	% testcase 10-12: large directed graphs

	% PARAMETERS: 
	switch testcase
		case 1
			graphname = "minnesota";
			type = "u";
			t = 1;
			alpha = 0.6;
			truesol_method = 'eig';			% uses eig function
			maxit1 = 150;					% polynomial
			maxit2 = 150;					% S&I Moret-Novati
			maxit3 = 120;					% S&I spectral
			maxitg = 50;					% rational EDS
		case 2
			graphname = "minnesota";
			type = "u";
			t = 10;
			alpha = 0.5;
			truesol_method = 'eig';			% uses eig function
			maxit1 = 150;					% polynomial
			maxit2 = 110;					% S&I Moret-Novati
			maxit3 = 70;					% S&I spectral
			maxitg = 45;					% rational EDS
		case 3
			graphname = "minnesota";
			type = "u";
			t = 1;
			alpha = 1;
			truesol_method = 'eig';			% uses eig function
			maxit1 = 25;					% polynomial
			maxit2 = 50;					% S&I Moret-Novati
			maxit3 = 130;					% S&I spectral
			maxitg = 70;					% rational EDS
		case 4
			graphname = "wiki-Vote";
			type = "d";
			t = 2;
			alpha = 0.25;
			truesol_method = 'eig';			% uses eig function
			maxit1 = 150;					% polynomial
			maxit2 = 150;					% S&I Moret-Novati
			maxit3 = 90;					% S&I spectral
			maxitg = 40;					% rational EDS
		case 5
			graphname = "wiki-Vote";
			type = "d";
			t = 0.3;
			alpha = 0.6;
			truesol_method = 'eig';			% uses eig function
			maxit1 = 150;					% polynomial
			maxit2 = 150;					% S&I Moret-Novati
			maxit3 = 90;					% S&I spectral
			maxitg = 40;					% rational EDS
		case 6
			graphname = "wiki-Vote";
			type = "d";
			t = 10;
			alpha = 0.75;
			truesol_method = 'eig';			% uses eig function
			maxit1 = 150;					% polynomial
			maxit2 = 40;					% S&I Moret-Novati
			maxit3 = 90;					% S&I spectral
			maxitg = 50;					% rational EDS
		case 7
			graphname = "Oregon-1";
			type = "u";
			t = 1;
			alpha = 0.6;
			truesol_method = 'projection';	% uses rational krylov with proj
			maxit1 = 150;					% polynomial
			maxit2 = 120;					% S&I Moret-Novati
			maxit3 = 150;					% S&I spectral
			maxitg = 50;					% rational EDS
		case 8
			graphname = "ca-HepPh";
			type = "u";
			t = 10;
			alpha = 0.6;
			truesol_method = 'projection';	% uses rational krylov with proj
			maxit1 = 150;					% polynomial
			maxit2 = 80;					% S&I Moret-Novati
			maxit3 = 110;					% S&I spectral
			maxitg = 50;					% rational EDS
		case 9
			graphname = "as-22july06";
			type = "u";
			t = 2;
			alpha = 0.4;
			truesol_method = 'projection';	% uses rational krylov with proj
			maxit1 = 150;					% polynomial
			maxit2 = 150;					% S&I Moret-Novati
			maxit3 = 150;					% S&I spectral
			maxitg = 50;					% rational EDS
		case 10
			graphname = "enron";
			type = "d";
			t = 10;
			alpha = 0.6;
			truesol_method = 'projection';	% uses rational krylov with proj
			maxit1 = 150;					% polynomial
			maxit2 = 50;					% S&I Moret-Novati
			maxit3 = 110;					% S&I spectral
			maxitg = 50;					% rational EDS
		case 11
			graphname = "p2p-Gnutella30";
			type = "d";
			t = 1;
			alpha = 0.6;
			truesol_method = 'projection';	% uses rational krylov with proj
			maxit1 = 130;					% polynomial
			maxit2 = 60;					% S&I Moret-Novati
			maxit3 = 50;					% S&I spectral
			maxitg = 30;					% rational EDS
		case 12
			graphname = "hvdc1";
			type = "d";
			t = 2;
			alpha = 0.4;
			truesol_method = 'projection';	% uses rational krylov with proj
			maxit1 = 200;					% polynomial
			maxit2 = 200;					% S&I Moret-Novati
			maxit3 = 180;					% S&I spectral
			maxitg = 50;					% rational EDS
		otherwise
			disp("unknown test case");
			return;
	end

	fv = @(B, w) fracexpmvec(B, alpha, t, w);		% exp(-t B^alpha)*w

	% UNDIRECTED TEST GRAPHS		nodes(LCC)		edges(LCC)
	% graphname = 'minnesota';		2640			6604			
	% graphname = 'Oregon-1';		11174			46818			
	% graphname = 'ca-HepPh';		11204			235238			
	% graphname = 'as-22july06'; 	22963			96872			

	% DIRECTED TEST GRAPHS				nodes(LCC)		edges(LCC)
	% graphname = 'wiki-Vote';			1300			39456
	% graphname = 'enron';				8271			146260		
	% graphname = 'p2p-Gnutella30';		8490			31706		
	% graphname = 'hvdc1';				24836			133620

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
	y1 = correction(y1, z);
	nerr1 = getErrorNorms(y_true, y1);

	% Shift-and-Invert: single repeated pole
	delta = t^(-2/alpha); % MORET-NOVATI CHOICE
	poles = -delta*ones(1, maxit2);
	[y2, t2] = runRatKrylov(L, v, poles, fv, "", maxit2);
	y2 = correction(y2, z);
	nerr2 = getErrorNorms(y_true, y2);

	delta = delta_spec; % SPECTRAL CHOICE
	poles = -delta*ones(1, maxit3);
	[y3, t3] = runRatKrylov(L, v, poles, fv, "", maxit3);
	y3 = correction(y3, z);
	nerr3 = getErrorNorms(y_true, y3);

	% General Krylov with asymptotically optimal poles (EDS)
	poles_eds = eds_cauchy_st(0.99*lambda_2, 1.01*lambda_n, maxitg+1);
	poles = poles_eds(1:maxitg);
	[yg, tg] = runRatKrylov(L, v, poles, fv, "", maxitg);
	yg = correction(yg, z);
	nerrg = getErrorNorms(y_true, yg);

	%%% ----- Rank-one shifted methods ----- %%%
	disp("Starting shifted methods");
	% Polynomial: all poles at infinity
	poles = Inf*ones(1, maxit1);
	[y1s, t1s] = runRatKrylov(L, v, poles, fv, "rk1shift", maxit1);
	nerr1s = getErrorNorms(y_true, y1s);

	% Shift-and-Invert: single repeated pole
	delta = t^(-2/alpha); % MORET-NOVATI CHOICE
	poles = -delta*ones(1, maxit2);
	[y2s, t2s] = runRatKrylov(L, v, poles, fv, "rk1shift", maxit2);
	nerr2s = getErrorNorms(y_true, y2s);

	delta = delta_spec; % SPECTRAL CHOICE
	poles = -delta*ones(1, maxit3);
	[y3s, t3s] = runRatKrylov(L, v, poles, fv, "rk1shift", maxit3);
	nerr3s = getErrorNorms(y_true, y3s);

	% General Krylov with asymptotically optimal poles (EDS)
	poles = poles_eds(1:maxitg);
	[ygs, tgs] = runRatKrylov(L, v, poles, fv, "rk1shift", maxitg);
	nerrgs = getErrorNorms(y_true, ygs);


	%%% ----- Projected methods ----- %%%
	disp("Starting projected methods");
	% Polynomial: all poles at infinity
	poles = Inf*ones(1, maxit1);
	[y1p, t1p] = runRatKrylov(L, v, poles, fv, "projection", maxit1);
	nerr1p = getErrorNorms(y_true, y1p);

	% Shift-and-Invert: single repeated pole
	delta = t^(-2/alpha); % MORET-NOVATI CHOICE
	poles = -delta*ones(1, maxit2);
	[y2p, t2p] = runRatKrylov(L, v, poles, fv, "projection", maxit2);
	nerr2p = getErrorNorms(y_true, y2p);

	delta = delta_spec;	% SPECTRAL CHOICE
	poles = -delta*ones(1, maxit3);
	[y3p, t3p] = runRatKrylov(L, v, poles, fv, "projection", maxit3);
	nerr3p = getErrorNorms(y_true, y3p);

	% General Krylov with asymptotically optimal poles (EDS)
	poles = poles_eds(1:maxitg);
	[ygp, tgp] = runRatKrylov(L, v, poles, fv, "projection", maxitg);
	nerrgp = getErrorNorms(y_true, ygp);

	%%% ----- Implicit projected methods ----- %%%
	disp("Starting implicit projected methods");
	% Compute f(L)v = f(L)w + f(L)z, where L z = 0 and w'z = 0
	w = v - z;
	fz = fv(0,1)*z; 
	% Polynomial: all poles at infinity
	poles = Inf*ones(1, maxit1);
	[y1ip, t1ip] = runRatKrylov(L, w, poles, fv, "", maxit1);
	y1ip = y1ip + fz;
	nerr1ip = getErrorNorms(y_true, y1ip);

	% Shift-and-Invert: single repeated pole
	delta = t^(-2/alpha); % MORET-NOVATI CHOICE
	poles = -delta*ones(1, maxit2);
	[y2ip, t2ip] = runRatKrylov(L, w, poles, fv, "", maxit2);
	y2ip = y2ip + fz;
	nerr2ip = getErrorNorms(y_true, y2ip);

	delta = delta_spec;	% SPECTRAL CHOICE
	poles = -delta*ones(1, maxit3);
	[y3ip, t3ip] = runRatKrylov(L, w, poles, fv, "", maxit3);
	y3ip = y3ip + fz;
	nerr3ip = getErrorNorms(y_true, y3ip);

	% General Krylov with asymptotically optimal poles (EDS)
	poles = poles_eds(1:maxitg);
	[ygip, tgip] = runRatKrylov(L, w, poles, fv, "", maxitg);
	ygip = ygip + fz;
	nerrgip = getErrorNorms(y_true, ygip);

	%%% ----------------------------------------------
	% Legend for naming: y<Krylov-type><desing-type>
	%	<Krylov-type>:	1 = polynomial
	%					2 = Moret-Novati S&I
	%					3 = Spectral S&I
	%					g = EDS Rational
	%	<desing-type>:	nothing = standard method
	%					s = rank-one shift
	%					p = projection
	%					ip = implicit projection
	%%% ----------------------------------------------

	% Display times:
	fprintf("Times for Krylov methods:\n");
	fprintf("\t\t\t Polynomial \t Moret-Novati SI \t Spectral SI \t EDS Rational\n");
	fprintf("Standard: \t\t %.3f \t\t %.3f \t\t\t %.3f \t\t %.3f \t\t\n", sum(t1), sum(t2), sum(t3), sum(tg));
	fprintf("Rank-one shift: \t %.3f \t\t %.3f \t\t\t %.3f \t\t %.3f \t\t\n", sum(t1s), sum(t2s), sum(t3s), sum(tgs));
	fprintf("Projection: \t\t %.3f \t\t %.3f \t\t\t %.3f \t\t %.3f \t\t\n", sum(t1p), sum(t2p), sum(t3p), sum(tgp));
	fprintf("Implicit projection: \t %.3f \t\t %.3f \t\t\t %.3f \t\t %.3f \t\t\n", sum(t1ip), sum(t2ip), sum(t3ip), sum(tgip));

	plottingScript;

end

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
