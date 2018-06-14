function fn = cell2node(x, fc, varargin)
xc = x(1:end-1) + 0.5*diff(x);

if nargin < 3
    method = 'linear';
else
    method = string(lower(varargin));
end

switch method
    case 'linear'
        fn = interp1(xc, fc, x, 'linear', 'extrap');
    case 'harmonic'
        fn = interp1(xc, 1./fc, x, 'linear', 'extrap');
        fn(isnan(fn)) = 0;
        fn(isinf(fn)) = 0;
        fn = 1./fn;
end
end