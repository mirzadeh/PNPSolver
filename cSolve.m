function cn = cSolve(x, dt, psi, c, cn, z, bc)
dx = diff(x);
f = getFlux(x, c, psi, z, 'linear');
f(1) = bc(1); f(end) = bc(end);
F = cn/dt - diff(f) ./ dx;

% note: we are solving c on cell centers
xc = x(1:end-1) + 0.5*dx;
A = 1/dt*speye(length(xc)) + matGen(x, 'cell');
cn = A \ F;

end
    