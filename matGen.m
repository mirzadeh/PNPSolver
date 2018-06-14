function A = matGen(x, varargin)
dx = diff(x);

if nargin < 2
    type = 'node';
else
    type = 'cell';
end

switch type
    case 'node'
        nx = length(x);
        dxl = [dx(1); dx];
        dxr = [dx; dx(end)];
        
        wl = -2./(dxl.*(dxl+dxr)); wl(1) = 0;
        wr = -2./(dxr.*(dxl+dxr)); wr(end) = 0;
        wc = -(wr+wl);
        
        % apply bc
        wc(1) = 1; wc(end) = 1;
        wr(1) = 0; wl(end) = 0;
        
    case 'cell'
        xc = x(1:end-1) + 0.5*dx;
        nx = length(xc);
        dxc = diff(xc);
        dxl = [dxc(1); dxc];
        dxr = [dxc; dxc(end)];
        
        wl = -1./(dxl.*dx); wl(1) = 0;
        wr = -1./(dxr.*dx); wr(end) = 0;
        wc = -(wr+wl);     
end

A = spdiags([wr wc wl], [-1 0 1], nx, nx)';
end