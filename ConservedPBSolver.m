function sol = ConservedPBSolver(x, lambda, v)
sol.grid = struct('x', x);
sol.options = struct('itmax', 100, 'tol', 1e-6, 'omega', 1);
sol.psi = zeros(size(x));
sol.err = [];

err = inf;
A = matGen(x);

sol.c0 = 1;
while err > sol.options.tol
    rhs = sol.c0/lambda^2*(sol.psi.*cosh(sol.psi) - sinh(sol.psi));
    rhs(1) = -v; rhs(end) = v;
    
    D = sol.c0/lambda^2 * sparse(diag(cosh(sol.psi)));
    D(1) = 0; D(end) = 0;
    psi_new = (A + D) \ rhs;
    
    err = max(abs(psi_new - sol.psi));
    sol.err = cat(1, sol.err, err);
    sol.psi = sol.psi + sol.options.omega*(psi_new - sol.psi);
    
    sol.c0 = 2/trapz(x, exp(sol.psi));
end
end