function [du] = grad(x,dx,u)
Nx = length(x);

du = zeros(Nx-1,1);

for i=1:Nx-1
%    du(i) = 0.5*((u(i+1)-u(i))/dx(i)+(u(i)-u(i-1))/dx(i-1));
    du(i) = (u(i+1)-u(i))/dx(i);
end

end
