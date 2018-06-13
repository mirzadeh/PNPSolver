function f = getFlux(x, c, psi, z, varargin)
if nargin < 5
    method = 'linear';
else
    method = lower(varargin{1});
end

e = -grad(psi, x, 'node');
f = z*e.*cell2node(x, c, method);
end