function ion = ion(z, c0, varargin)
    if nargin < 3
        d = 1;
    else
        d = varargin{1};
    end
    
    ion = struct('z', z, 'c0', c0, 'd', d);
end