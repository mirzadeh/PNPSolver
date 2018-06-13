function sol = PNPSolver(x, lambda, tf)

dx = diff(x);
xc = x(1:end-1) + 0.5*dx;
sol.grid = struct('x', x, 'xc', xc, 'dx', dx);
sol.options = struct('lambda', lambda, 'tf', tf, 'iter_max', 1, 'tol', 1e-6, ...
    'dt', 0.05*lambda^2);

nx = length(x);
sol.cp  = 0.5*ones(nx-1,1);
sol.cm  = 0.5*ones(nx-1,1);
sol.psi = zeros(nx,1);

sol.q = [];
sol.t = [];

dt  = sol.options.dt;
t   = 0;
tc  = 0;

cp_n  = sol.cp;
cm_n  = sol.cm;
psi_n = sol.psi;

while(t < tf)    
    iter = 1;
    err = 1;
    
    fprintf(' ---------------- Iteration ---------------- \n');
    fprintf(' t = %f \t time-step = %d\n\n', t, tc);
    cp_tmp  = cp_n;
    cm_tmp  = cm_n;
    psi_tmp = psi_n;
    while (iter <= sol.options.iter_max && err > sol.options.tol)
        cp_tmp  = cSolve(x, dt, psi_tmp, cp_tmp, cp_n,  1, [0 0]);
        cm_tmp  = cSolve(x, dt, psi_tmp, cm_tmp, cm_n, -1, [0 0]);        
        psi_tmp = pSolve(x, cp_tmp, cm_tmp, 1/lambda, [-0.1 0.1]);

        d_cp  = max(abs(cp_n  - cp_tmp)); 
        d_cm  = max(abs(cm_n  - cm_tmp)); 
        d_psi = max(abs(psi_n - psi_tmp));
        
        iter = iter + 1;
        err = max(d_cp,max(d_cm,d_psi));
        
        fprintf(' iter = %2d \t err = %e\n',iter,err);
    end
    cp_n  = cp_tmp;
    cm_n  = cm_tmp;
    psi_n = psi_tmp;
    fprintf(' ---------------- ********* ---------------- \n');

    sol.cp  = cat(2, sol.cp, cp_n);
    sol.cm  = cat(2, sol.cm, cm_n);
    sol.psi = cat(2, sol.psi, psi_n);

    sol.t   = cat(1, sol.t, t);
    sol.q   = cat(1, sol.q, integrate(x, cp_n - cm_n));
    
    tc = tc + 1; t = t + dt;    
end
end