function [cpb,cmb,psib] = PB(x,cp,cm,psi)

Nx  = length(x);
[y] = MAC(x);

psib = cp;

for i=1:Nx-1
    psib(i) = inter(x(i),x(i+1),psi(i),psi(i+1),y(i));
end

if (mod(Nx,2)==0)
    kcp = cp(Nx/2);
    kcm = cm(Nx/2);
else
    im  = (Nx-1)/2;
    ip  = (Nx+1)/2;
    kcp = inter(x(im),x(ip),cp(im),cp(ip),0.5*(x(ip)+x(im)));
    kcm = inter(x(im),x(ip),cm(im),cm(ip),0.5*(x(ip)+x(im)));
end

cpb = kcp*exp(-psib);
cmb = kcm*exp( psib);
end