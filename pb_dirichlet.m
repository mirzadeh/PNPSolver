function [cp, cm, psi] = pb_dirichlet(v, k, x)
    [cp, cm, psi] = solve_nonlinear(x, v, k);
end

function [cp, cm, psi] = solve_nonlinear(x, v, k)
    dx = diff(x);
    n = length(x);
    psi_old = zeros(n,1);
    psi_new = psi_old;
    
    B = zeros(n, 3);
    
    err = 1;
    tol = 1e-6;
    
    B(1,2) = 1;
    for i=2:n-1
        B(i-1,1) = -2./(dx(i-1)*(dx(i)+dx(i-1)));
        B(i+1,3) = -2./(dx(i)  *(dx(i)+dx(i-1)));
        B(i,  2) = -(B(i-1,1) + B(i+1,3));
    end
    B(n,2) = 1;

    c0 = 1;
    while err > tol
        c0 = 2/trapz(x,exp(psi_old));
        diags = B;
        diags(2:n-1,2) = B(2:n-1,2) + c0*k^2*cosh(psi_old(2:n-1));
                
        A = spdiags(diags,[-1,0,1],n,n);                
        rhs = c0*k^2*(psi_old.*cosh(psi_old) - sinh(psi_old));
        rhs(1) = -v; rhs(n) = v;

        psi_new = A \ rhs;        
        
        err = max(abs(psi_new - psi_old));        
        psi_old = psi_old + 0.5*(psi_new-psi_old);
    end
    disp(err)
    
    psi = psi_new;            
    cp = 0.5*c0*exp(-psi);
    cm = 0.5*c0*exp( psi);
end
