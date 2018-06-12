function [B] = vecGen(x,dx,bcv,f,Type)
Nx = length(x);
B  = zeros(Nx,1);

switch Type
    case 'Dirichlet'
        B(1) = bcv(1);
        for i=2:Nx-1
            B(i) = f(i);
        end
        B(Nx) = bcv(2); 
    case 'Neumann'
        B(1) = f(1) - 2*bcv(1)/dx(1);
        for i=2:Nx-1
            B(i) = f(i);
        end
        B(Nx) = f(Nx) + 2*bcv(2)/dx(Nx-1);
    case 'NP'
        B(1) = f(1);
        for i=2:Nx-1
            B(i) = f(i);
        end
        B(Nx) = f(Nx);

end
end