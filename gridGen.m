function [x, dx] = gridGen(xmin, xmax, Nx, varargin)
if (nargin < 4 )
    x = linspace(xmin, xmax, Nx);
else
    tau = varargin{1};
    x = linspace(-1, 1, Nx);
    A = 1/(erf(tau) - erf(-tau));
    x = xmin + (xmax - xmin)*A*(erf(tau*x) - erf(-tau));
end

dx = diff(x);
end