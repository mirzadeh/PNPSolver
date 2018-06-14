function fc = node2cell(x, fn, varargin)
xc = x(1:end-1) + 0.5*diff(x);

if nargin < 3
    method = 'linear';
else
    method = string(lower(varargin));
end

switch method
    case 'linear'
        fc = interp1(x, fn, xc, 'linear', 'extrap');
    case 'harmonic'
        fc = interp1(x, 1./fn, xc, 'linear', 'extrap');
        fc(isnan(fc)) = 0;
        fc(isinf(fc)) = 0;
        fc = 1./fc;
end
end