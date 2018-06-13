function x = gridGen(xmin, xmax, nx, varargin)
if (nargin < 4 )
    x = linspace(xmin, xmax, nx)';
else
    tau = varargin{1};
    x = linspace(-1, 1, nx)';
    A = 1/(erf(tau) - erf(-tau));
    x = xmin + (xmax - xmin)*A*(erf(tau*x) - erf(-tau));
end
end