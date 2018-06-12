function [f] = getFlux(x,dx,c,psi,z,d)
Nx = length(x);
f  = zeros(Nx,1);
elec = -grad(x,dx,psi);
[y,dy] = MAC(x);

for i=2:Nx-1
    f(i) = z * 0.5*d*(c(i)*elec(i)+c(i-1)*elec(i-1));
    %f(i) = f(i) + d*(c(i)-c(i-1))/dy(i-1);
end

end