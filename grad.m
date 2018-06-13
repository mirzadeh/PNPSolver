function du = grad(u, x, varargin)
dx = diff(x);
if nargin < 3
    type = 'node';
else
    type = lower(varargin{1});
end

switch type
    case 'node'
        du = cell2node(x, diff(u) ./ dx);
    case 'cell'
        xc = x(1:end-1) + 0.5*dx;
        du = interp1(x(2:end-1), diff(u) ./ diff(xc), x, 'linear', 'extrap');
end
