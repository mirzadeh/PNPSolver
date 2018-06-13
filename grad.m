function du = grad(u, x, dx, varargin)
if nargin < 3
    type = 'node';
else
    type = lower(varargin{1});
end

switch type
    case 'node'
        du = cell2node(diff(u) ./ dx);
    case 'cell'
        xc = x(1:end-1) + 0.5*dx;
        du = diff(u) ./ diff(xc);
        dul = inter(x(2), x(3), du(1), du(2), x(1));
        dur = inter(x(end-1), x(end-2), du(end), du(end-1), x(end));
        du = [dul, du, dur];
end
