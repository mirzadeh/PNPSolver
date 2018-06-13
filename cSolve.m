function cn = cSolve(x, dt, psi, c, cn, z, bc)
dx = diff(x);
f = getFlux(x, c, psi, z, 'linear');
F = cn/dt - diff(f) ./ dx;

% apply flux boundary conditions
F(1) = F(1) + bc(1)/dx(1);
F(end) = F(end) - bc(end)/dx(end);

% note: we are solving c on cell centers
xc = x(1:end-1) + 0.5*diff(x);
A = 1/dt*speye(length(xc)) + matGen(xc, 'neumann');
cn = A \ F;

end
    