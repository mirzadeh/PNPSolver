function [sum] = nint(x,dx,inter,v)
Nx = length(x);
sum = 0;
for i=1:Nx-1
    if ( x(i) < inter(1) ) continue; end
    if ( x(i) >= inter(2) ) break; end
    
    sum = sum + v(i) * dx(i);
end
end