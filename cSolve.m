function cn = cSolve(x, dt, psi, c, z, d, bc)
dx = diff(x);
f = getFlux(dx, c, psi, d, z, 'linear');
F = c/dt - diff(f) ./ dx;

% apply flux boundary conditions
F(1) = F(1) + bc(1)/dx(1);
F(end) = F(end) - bc(end)/dx(end);

% note: we are solving c on cell centers
xc = x + 0.5*diff(x);
A = matGen(xc, 'neumann');
cn = A \ F;

end
    