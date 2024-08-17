function [y, tottime] = runRatKrylov(A, b, poles, fv, option, maxit)
	% [y, tottime] = runRatKrylov(A, b, poles, fv, option, maxit)
	% Computes the first m = min(length(poles), maxit) iterations of a rational Krylov method
	% fv(M, x) computes the product f(M) x for a matrix M and vector x
	% If A is the transpose of a graph Laplacian, the following options 
	% are available:
	% option = "rk1shift" perform rank one shift 
	% option = "projection" perform projection on ones(n,1)^\perp
	% anything else simply runs the Krylov method
	%
	% Output:
	% y(:, k) contains the approximation to f(A)*b after k iterations
	% tottime contains the total runtime for the algorithm 
	% (excluding preprocessing for the projection or rank-one shift)

	if nargin < 5
		option = "";
	end
	if nargin < 6
		maxit = length(poles);
	end
	n = size(A, 1);
	maxit = min(maxit, length(poles));
	y = zeros(n, maxit);
	if option == "rk1shift"		% using rank-one shift method
		if (max(max(abs(A-A'))) == 0)	% symmetric
			[x, z] = prepareShift(A, "u");
		else
			[x, z] = prepareShift(A);
		end
		tic;
		% x = ones(n,1) and sum(z) = 1, A z = 0
		% A_rk1shift = A + theta*z*x';	% exact shifted Laplacian 
		% f(A)b = f(A + zx')b - [f(theta)-f(0)]z x'b
		% backshift is such that f(A) b = f(A_rk1shift) b - backshift
		theta = 1;
		backshift = (fv(theta,1)-fv(0,1))*z*(x'*b);

		% run polynomial Krylov
		A_rk1shift.multiply = @(rho, eta, w) ...
			rho*A*w - eta*w + rho*theta*(x'*w)*z;
		A_rk1shift.isreal = isreal(A);
		if (min(poles) == max(poles) && poles(1) ~= Inf)
			% Prepare Shift-and-Invert Krylov:
			% Compute LU factorization of A - xi I only once
			xi = poles(1);
			perm = amd(A);
			[l, u] = lu(A(perm,perm) - xi*speye(n));
			A_rk1shift.solve = @(rho, eta, w) ...
				(1/rho)*( LUsolver(l,u,w - (x'*w)*z, perm) + ...
				(x'*w)/(-xi+theta)*z );
			% A_rk1shift.solve only works if the poles are all 
			% equal to xi and xi = eta/rho
		elseif min(poles) == Inf
			% prepare dummy solver (otherwise we may get Inf or NaN):
			A_rk1shift.solve = @(rho, eta, w) w;
		else
			% prepare standard shifted solver:
			% ( rho*(A + theta*z*x') - eta*speye(n) )\w;
			A_rk1shift.solve = @(rho, eta, w) ...
				1/rho*( (A - eta/rho*speye(n))\w + ...
				theta/((eta/rho)*(theta-eta/rho))*(x'*w)*z );
		end
		% run shifted Krylov algorithm
		[V, ~, ~, ~] = rat_krylov(A_rk1shift, b, poles(1:maxit));

		Am = V(:,1:maxit)'*A*V(:,1:maxit) + theta*(V(:,1:maxit)'*z)*(x'*V(:,1:maxit));
		for k = 1:maxit
			y(:,k) = V(:,1:k)*fv(Am(1:k,1:k), V(:,1:k)'*b);
		end

		y = y - backshift;
	elseif option == "projection"		% using projection method
		%%% Preparing projection: we compute f(L)b = f(L) w + f(L) z 
		% note that f(L)*z = f(0)*z
		if (max(max(abs(A-A'))) == 0)	% symmetric
			[x, z] = prepareShift(A, "u");
		else
			[x, z] = prepareShift(A);
		end
		tic;
		s = (-1-1/sqrt(n))/(n-1);
		w = b - z;		% w has zero sum
		
		% vv \in \R^{n-1}
		Qmult = @(vv) [1/sqrt(n)*(ones(n-1, 1)'*vv); vv + s*ones(n-1,1)*(ones(n-1, 1)'*vv)];
		% u \in \R^n
		QTmult = @(u) u(2:end,:) + ones(n-1,1)*( 1/sqrt(n)*u(1,:) + s*sum(u(2:end,:)) );
		QTw = QTmult(w); % starting vector for projected method
		QTAQ.multiply = @(rho, eta, xx) ...
			rho*QTmult(A*Qmult(xx)) - eta*xx;
		QTAQ.isreal = isreal(A);
		
		if (min(poles) == max(poles) && poles(1) ~= Inf)
			% prepare Shift-and-Invert Krylov solver:
			xi = poles(1);
			perm = amd(A);
			[l, u] = lu(A(perm, perm) - xi*speye(n));
			QTAQ.solve = @(rho, eta, xx) ...
				1/rho*QTmult(LUsolver(l, u, Qmult(xx), perm));
			else
			% prepare standard Krylov solver: 
			QTAQ.solve = @(nu, mu, xx) ...
				QTmult( (nu*A - mu*speye(n))\(Qmult(xx)) );
		end
		[V, ~, ~, ~] = rat_krylov(QTAQ, QTw, poles(1:maxit));

		QVm = Qmult(V(:, 1:maxit));

		projAm = QVm'*A*QVm;
		for k = 1:maxit
			y(:,k) = Qmult(V(:,1:k)*fv(projAm(1:k,1:k),V(:,1:k)'*QTw));
		end

		y = y + fv(0,1)*z;
	else 		% using "original" Krylov method
		tic;
		Astruct.multiply = @(rho, eta, w) rho*A*w - eta*w;
		Astruct.isreal = isreal(A);

		if (min(poles) == max(poles) && poles(1) ~= Inf)
			% prepare Shift-and-Invert Krylov:
			% Compute LU factorization of A - xi I only once
			xi = poles(1);
			perm = amd(A);
			[l, u] = lu(A(perm, perm) - xi*speye(n));
			Astruct.solve = @(rho, eta, w) 1/rho*LUsolver(l,u,w, perm);
		else
			% prepare standard solver
			Astruct.solve = @(rho, eta, w) (rho*A - eta*speye(n))\w;
		end
		% run rational Krylov
		[V, ~, ~, ~] = rat_krylov(Astruct, b, poles(1:maxit));

		Am = V(:,1:maxit)'*A*V(:,1:maxit);
		for k = 1:maxit
			y(:,k) = V(:,1:k)*fv(Am(1:k,1:k), V(:,1:k)'*b);
		end

	end
	tottime = toc;
end



