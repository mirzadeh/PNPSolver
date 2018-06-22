function f  = laplace(x, u, varargin)
    dx = diff(x);
    xc = x(1:end-1) + 0.5*dx;
    if nargin < 3
        type = 'node';
    else
        type = lower(varargin{1});
    end
    
    f = zeros(size(u));
    switch type
        case 'node'
            f(2:end-1) = diff(diff(u) ./ dx) ./ diff(xc);
            f(1) = interp1(x(2:end-1), f(2:end-1), x(1), 'linear', 'extrap');
            f(end) = interp1(x(2:end-1), f(2:end-1), x(end), 'linear', 'extrap');
            
        case 'cell'            
            f = diff(diff(u) ./ diff(xc)) ./ dx(2:end-1);
            f(1) = interp1(xc(2:end-1), f(2:end-1), xc(1), 'linear', 'extrap');
            f(end) = interp1(xc(2:end-1), f(2:end-1), xc(end), 'linear', 'extrap');
    end
    
end