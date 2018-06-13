function sol = PNP_discharge(init, delta, beta, tf)

sol.runData.init      = init;
sol.runData.delta     = delta;
sol.runData.beta      = beta;
sol.runData.tf        = tf;
sol.runData.dt_factor = 0.5;
sol.runData.iter_max  = 1;
sol.runData.tol_iter  = 1e-6;
sol.runData.xmin      = -1;
sol.runData.xmax      =  1;
sol.runData.Nx        = 200;
sol.runData.taux      = 2;

dt = sol.runData.dt_factor*delta^2;
dt1 = dt;
iter_max = sol.runData.iter_max;
tol_iter = sol.runData.tol_iter;

[x,dx] = gridGen(sol.runData.xmin, sol.runData.xmax, sol.runData.Nx, sol.runData.taux);
[y,~]  = MAC(x);

psi = init.psi;
cp = init.cp;
cm = init.cm;

sol.cp = cp;
sol.cm = cm;
sol.psi = psi;

sol.t = 0;
sol.x = x;
sol.y = y;
sol.Sigma = 0;

t   = 0;
tc  = 0;
dt = dt1;

bcv = [0 0];
sol.q = nint(x, dx, [-1,0], cp-cm);

sol.frames.q(1000) = struct('cdata',[],'colormap',[]);

while(t<sol.runData.tf)
    
    iter = 1;
    err_iter = 1;
    fprintf(' ---------------- Iteration ---------------- \n');
    fprintf(' t = %f \t time-step = %d\n\n', t, tc);
    cp_tmp  = cp;
    cm_tmp  = cm;
    while (iter<=iter_max && err_iter>tol_iter)
        [cp_n]  = cSolve(x,dx,dt,psi,cp, 1, 1);
        [cm_n]  = cSolve(x,dx,dt,psi,cm,-1, beta);
        [psi_n] = pSolve(x,dx,cp_n,cm_n,1/delta,bcv);

        d_cp  = max(abs(cp_n-cp_tmp)); cp_tmp = cp_n;  
        d_cm  = max(abs(cm_n-cm_tmp)); cm_tmp = cm_n; 
        d_psi = max(abs(psi -psi_n));     psi = psi_n;
        
        iter = iter + 1;
        err_iter = max(d_cp,max(d_cm,d_psi));
        
        fprintf(' iter = %2d \t err = %e\n',iter,err_iter);
    end
    fprintf(' ---------------- ********* ---------------- \n');

    dt = min(1.01*dt,dt1);
    
    cp = cp_n;
    cm = cm_n;
            
    tc = tc + 1; t = t + dt;
    
    q = nint(x, dx, [-1,0], cp-cm);
    
    sol.cp  = cat(2, sol.cp, cp);
    sol.cm  = cat(2, sol.cm, cm);
    sol.psi = cat(2, sol.psi, psi);

    sol.t   = cat(2, sol.t, t);
    sol.q   = cat(2, sol.q, q);
    
    sol.Sigma = cat(2, sol.Sigma, -nint(x,dx, [sol.runData.xmin 0], (cp-cm))*delta);

end
end