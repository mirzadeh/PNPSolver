function sol = PBSolver(x, lambda, v)
sol.grid = struct('x', x);
sol.options = struct('itmax', 100, 'tol', 1e-6, 'omega', 1);
sol.psi = zeros(size(x));
sol.err = [];

err = inf;
A = matGen(x);
while err > sol.options.tol
    rhs = 1/lambda^2*(sol.psi.*cosh(sol.psi) - sinh(sol.psi));
    rhs(1) = -v; rhs(end) = v;
    
    D = 1/lambda^2 * sparse(diag(cosh(sol.psi)));
    D(1) = 0; D(end) = 0;
    psi_new = (A + D) \ rhs;
    
    err = max(abs(psi_new - sol.psi));
    sol.err = cat(1, sol.err, err);
    sol.psi = sol.psi + sol.options.omega*(psi_new - sol.psi);
end
end