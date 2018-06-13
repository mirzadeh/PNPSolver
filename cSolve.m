function cn = cSolve(dx, dt, psi, c, z, d, bc)
f = getFlux(dx, c, psi, z, d, 'linear');
F = c/dt - diff(f) ./ dx;

% apply flux boundary conditions
F(1) = F(1) + bc(1)/dx(1);
F(end) = F(end) - bc(end)/dx(end);

A = matGen(y,dy,0,0,dt,'NP', d);
cn = A\F;

end
    