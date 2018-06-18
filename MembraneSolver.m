function sol = MembraneSolver(x, lambda, v, G, lambda_m, tf)
dx.l = diff(x.l);
xc.l = x.l(1:end-1) + 0.5*dx.l;

dx.r = diff(x.r);
xc.r = x.r(1:end-1) + 0.5*dx.r;

sol.grid = struct('xl', x.l, 'xr', x.r, 'xcl', xc.l, 'xcr', xc.r, 'dxl', dx.l, 'dxr', dx.r);
sol.options = struct('lambda', lambda, 'tf', tf, 'v', v, 'G', G, ...
    'lambda_m', lambda_m, 'iter_max', 5, 'tol', 1e-6, 'dtmax', 0.5*lambda^2);

nl = length(x.l);
cp_n.l  = 0.5*ones(nl-1,1);
cm_n.l  = 0.5*ones(nl-1,1);
psi_n.l = linspace(-v, v, nl)';

nr = length(x.r);
cp_n.r  = 0.5*ones(nr-1,1);
cm_n.r  = 0.5*ones(nr-1,1);
psi_n.r = linspace(-v, v, nr)';

sol.t = 0;
sol.cp.l = cp_n.l;
sol.cm.l = cm_n.l;
sol.psi.l = psi_n.l;
sol.cp.r = cp_n.r;
sol.cm.r = cm_n.r;
sol.psi.r = psi_n.r;


t   = 0;
tc  = 0;
dt = min(0.5*min([sol.grid.dxl; sol.grid.dxr]), sol.options.dtmax);

while(t < tf)    
    iter = 1;
    err = 1;
    
    fprintf(' ---------------- Iteration ---------------- \n');
    fprintf(' t = %f \t time-step = %d\n\n', t, tc);
    cp_tmp  = cp_n;
    cm_tmp  = cm_n;
    psi_tmp = psi_n;
    while (iter <= sol.options.iter_max && err > sol.options.tol)        
        cp_new  = cSolve(x, dt, psi_tmp, cp_tmp, cp_n,  1, G);
        cm_new  = cSolve(x, dt, psi_tmp, cm_tmp, cm_n, -1, G);                
        psi_new = pSolve(x, cp_tmp, cm_tmp, 1/lambda, lambda_m, [-v v]);

        d_cp.l  = integrate(x.l, (cp_new.l  - cp_tmp.l).^2); 
        d_cm.l  = integrate(x.l, (cm_new.l  - cm_tmp.l).^2); 
        d_psi.l = integrate(x.l, (psi_new.l - psi_tmp.l).^2, 'node');
        
        d_cp.r  = integrate(x.r, (cp_new.r  - cp_tmp.r).^2); 
        d_cm.r  = integrate(x.r, (cm_new.r  - cm_tmp.r).^2); 
        d_psi.r = integrate(x.r, (psi_new.r - psi_tmp.r).^2, 'node');
        
        iter = iter + 1;
        err = max([d_cp.l d_cp.r d_cm.l d_cm.r d_psi.l d_psi.r]);
        
        cp_tmp = cp_new;
        cm_tmp = cm_new;
        psi_tmp = psi_new;
        fprintf(' iter = %2d \t err = %e\n', iter, err);
    end
    cp_n  = cp_tmp;
    cm_n  = cm_tmp;
    psi_n = psi_tmp;
    fprintf(' ---------------- ********* ---------------- \n');

    sol.cp.l  = cat(2, sol.cp.l, cp_n.l);
    sol.cm.l  = cat(2, sol.cm.l, cm_n.l);
    sol.psi.l = cat(2, sol.psi.l, psi_n.l);
    
    sol.cp.r  = cat(2, sol.cp.r, cp_n.r);
    sol.cm.r  = cat(2, sol.cm.r, cm_n.r);
    sol.psi.r = cat(2, sol.psi.r, psi_n.r);
        
    dt = min(1.02*dt, sol.options.dtmax);
    tc = tc + 1; t = t + dt;    
    sol.t = cat(1, sol.t, t);
end
end

function pn = pSolve(x, cp, cm, kappa, lambda_m, bc)
% since we are solving the Poisson equation using jump, we need to construct
% individual blocks separately
%
% A = [Al Alr; Arl Ar]
% left  equation: p_r - p_l = lambda_m * (p_l - p_{l-1})/dxl
% right equation: p_r - p_l = lambda_m * (p_{r+1} - p_r)/dxr
nl = length(x.l);
nr = length(x.r);
dxl = x.l(end) - x.l(end-1);
dxr = x.r(2) - x.r(1);

Al = matGen(x.l);
Ar = matGen(x.r);

Al(nl, nl-1) = lambda_m/dxl;
Al(nl, nl) = -1 - lambda_m/dxl;
Ar(1, 1) = 1 + lambda_m/dxr;
Ar(1, 2) = -lambda_m/dxr;

% include the off-diagonal terms
Al_off = sparse(nl, nr);
Ar_off = sparse(nr, nl);
Al_off(nl, 1) =  1;
Ar_off(1, nr) = -1;

A = [[Al Al_off]; [Ar_off Ar]]; 

% interpolate charge density on nodes
f.l = cell2node(x.l, kappa^2*(cp.l - cm.l));
f.r = cell2node(x.r, kappa^2*(cp.r - cm.r));

% adjust boundary conditions
f.l(1) = bc(1);
f.r(end) = bc(end);

sol = A \ [f.l; f.r];

pn.l = sol(1:nl);
pn.r = sol(nl+1:end);
end

function cn = cSolve(x, dt, psi, c, cn, z, G)
    % compute jmem -- only selective to cations.
    switch z
        case 1
            jmem = -G*(log(c.r(1)/c.l(end)) + (psi.r(1) - psi.l(end)));
        case -1
            jmem = 0;
    end

    % solve on either side using jmem as boundary condition
    cn.l = solveOneSide(x.l, dt, psi.l, c.l, cn.l, z, [0 jmem]);
    cn.r = solveOneSide(x.r, dt, psi.r, c.r, cn.r, z, [jmem 0]);

    function cn = solveOneSide(x, dt, psi, c, cn, z, bc)
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
end