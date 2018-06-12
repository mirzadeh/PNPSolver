function [x,dx] = gridGen(xmin, xmax, Nx, tau)
xtemp = linspace(xmin,xmax,Nx);

if (tau < 0 )
    x = xtemp;
else
ztemp = linspace(-1,1,Nx);
A = 2/(erf(tau) - erf(-tau));
B = -1 - A*erf(-tau);
zeta = A*erf(tau*ztemp) + B;
zeta  = zeta*(xmax-xmin)/2;
x  = zeta + xmin - zeta(1);

end
dx = zeros(Nx-1,1);
for i=1:Nx-1
    dx(i) = x(i+1) - x(i);
end
end