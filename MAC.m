function [z,dz] = MAC(x)

Nx = length(x);
z  = zeros(Nx-1,1);
dz = zeros(Nx-2,1);
for i=1:Nx-1
    z(i) = 0.5*(x(i)+x(i+1));
end
for i=1:Nx-2
    dz(i) = z(i+1)-z(i);
end
end