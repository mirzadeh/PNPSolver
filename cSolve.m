function [cn] = cSolve(x,dx,dt,psi,c,z,d)
Nx = length(x);
[y,dy] = MAC(x);

f = getFlux(x,dx,c,psi,z,d); F=y;

F(1) = c(1)/dt - (f(2)-f(1))*2/(dy(1)+dy(1));
for i=2:Nx-2
    F(i) = c(i)/dt - (f(i+1)-f(i))*2/(dy(i)+dy(i-1));
end
F(Nx-1) = c(Nx-1)/dt - (f(Nx)-f(Nx-1))*2/(dy(Nx-2)+dy(Nx-2));

A = matGen(y,dy,0,0,dt,'NP', d);
cn = A\F;

end
    