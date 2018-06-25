function ion = ion(z, c0, G, d)
    if nargin < 3
        G = 0;
    end
    if nargin < 4
        d = 1;
    end
    
    ion = struct('z', z, 'c0', c0, 'G', G, 'd', d);
end