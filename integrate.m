function res = integrate(x, f, varargin)
if nargin < 3
    type = 'cell';
else
    type = lower(varargin{1});
end

dx = diff(x);
switch type
    case 'cell'
        res = sum(dx.*f);
    case 'node'
        res = trapz(x, f);
end
end