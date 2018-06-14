function sol = PNPSolver(x, lambda, v, tf)
dx = diff(x);
xc = x(1:end-1) + 0.5*dx;
sol.grid = struct('x', x, 'xc', xc, 'dx', dx);
sol.options = struct('lambda', lambda, 'tf', tf, 'v', v, ...
    'iter_max', 5, 'tol', 1e-6, 'dtmax', 0.5*lambda^2);

nx = length(x);
cp_n  = 0.5*ones(nx-1,1);
cm_n  = 0.5*ones(nx-1,1);
psi_n = linspace(-v, v, nx)';

sol.q = 0;
sol.t = 0;
sol.cp = cp_n;
sol.cm = cm_n;
sol.psi = psi_n;

t   = 0;
tc  = 0;
dt = min(0.5*min(dx), sol.options.dtmax);

while(t < tf)    
    iter = 1;
    err = 1;
    
    fprintf(' ---------------- Iteration ---------------- \n');
    fprintf(' t = %f \t time-step = %d\n\n', t, tc);
    cp_tmp  = cp_n;
    cm_tmp  = cm_n;
    psi_tmp = psi_n;
    while (iter <= sol.options.iter_max && err > sol.options.tol)        
        cp_new  = cSolve(x, dt, psi_tmp, cp_tmp, cp_n,  1, [0 0]);
        cm_new  = cSolve(x, dt, psi_tmp, cm_tmp, cm_n, -1, [0 0]);                
        psi_new = pSolve(x, cp_tmp, cm_tmp, 1/lambda, [-v v]);

        d_cp  = integrate(x, (cp_new  - cp_tmp).^2); 
        d_cm  = integrate(x, (cm_new  - cm_tmp).^2); 
        d_psi = integrate(x, (psi_new - psi_tmp).^2, 'node');
        
        iter = iter + 1;
        err = max(d_cp, max(d_cm, d_psi));
        
        cp_tmp = cp_new;
        cm_tmp = cm_new;
        psi_tmp = psi_new;
        fprintf(' iter = %2d \t err = %e\n',iter,err);
    end
    cp_n  = cp_tmp;
    cm_n  = cm_tmp;
    psi_n = psi_tmp;
    fprintf(' ---------------- ********* ---------------- \n');

    sol.cp  = cat(2, sol.cp, cp_n);
    sol.cm  = cat(2, sol.cm, cm_n);
    sol.psi = cat(2, sol.psi, psi_n);
    
    mid = floor(nx/2);
    sol.q   = cat(1, sol.q, integrate(x(1:mid+1), cp_n(1:mid) - cm_n(1:mid)));
    
    dt = min(1.02*dt, sol.options.dtmax);
    tc = tc + 1; t = t + dt;    
    sol.t = cat(1, sol.t, t);
end
end

function pn = pSolve(x, cp, cm, kappa, bc)
% interpolate charge density on nodes
f = cell2node(x, kappa^2*(cp-cm));

% adjust boundary conditions
f(1) = bc(1);
f(end) = bc(end);

pn = matGen(x) \ f;
end

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

function f = getFlux(x, c, psi, z, varargin)
if nargin < 5
    method = 'linear';
else
    method = lower(varargin{1});
end

e = -grad(x, psi);
f = z*e.*cell2node(x, c, method);
end