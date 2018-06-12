function [pn] = pSolve(x,dx,cp,cm,kappa,bcv)
Nx = length(x);
[y,dy]=MAC(x);

F = zeros(Nx,1);
f = kappa*kappa*(cp-cm);


for i=2:Nx-1
    F(i) = inter(y(i-1),y(i),f(i-1),f(i),x(i));
end


A = matGen(x,dx,0,0,0,'Dirichlet', 1);
G = vecGen(x,dx,bcv,F,'Dirichlet');

pn = A\G;
end