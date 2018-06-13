function f = getFlux(x, c, psi, d, z, varargin)
dx = diff(x);

if nargin < 6
    method = 'linear';
else
    method = lower(varargin{1});
end

e = -grad(psi, x, dx, 'node');
f = z*d*e.*cell2node(c, method);
end