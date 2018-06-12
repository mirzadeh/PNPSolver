function [A] = matGen(x,dx,E,z,dt,Type, d)
Nx = length (x);
A  = sparse (Nx,Nx);
switch Type
    case 'Dirichlet'
        A(1,1) = 1;
        for i=2:Nx-1
            A(i,i-1) = -2./(dx(i-1)*(dx(i)+dx(i-1)));
            A(i,i)   =  2./(dx(i-1)* dx(i));
            A(i,i+1) = -2./(dx(i)  *(dx(i)+dx(i-1)));
        end
        A(Nx,Nx) = 1.;
    case 'Neumann'
        A(1,1) =  2.*d/dx(1)^2;
        A(1,2) = -2.*d/dx(1)^2;
        for i=2:Nx-1
            A(i,i-1) = -2.*d/(dx(i-1)*(dx(i)+dx(i-1)));
            A(i,i)   =  2.*d/(dx(i-1)* dx(i));
            A(i,i+1) = -2.*d/(dx(i)  *(dx(i)+dx(i-1)));
        end
        A(Nx,Nx-1) = -2.*d/dx(Nx-1)^2;
        A(Nx,Nx)   =  2.*d/dx(Nx-1)^2;
        
    case 'NP'
        A(1,1) =  1./dt + 1.*d/(dx(1)^2);
        A(1,2) = -1.*d/dx(1)^2;
        for i=2:Nx-1
            A(i,i-1) = -2.*d/(dx(i-1)*(dx(i)+dx(i-1)));
            A(i,i)   =  1./dt + 2.*d/(dx(i-1)* dx(i));
            A(i,i+1) = -2.*d/(dx(i)  *(dx(i)+dx(i-1)));
        end
        A(Nx,Nx-1) = -1.*d/dx(Nx-1)^2;
        A(Nx,Nx)   =  1./dt + 1.*d/dx(Nx-1)^2;
        
end
end