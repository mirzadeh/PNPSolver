function A = matGen(x, varargin)
dx = diff(x);
nx = length(x);

if nargin < 2
    bc = ["neumann", "neumann"];
else
    bc = string(lower(varargin{1}));
end

dxl = [dx(1) dx];
dxr = [dx dx(end)];
wl = -2./(dxl.*(dxl+dxr)); wl(1) = 0;
wr = -2./(dxr.*(dxl+dxr)); wr(end) = 0;
wc = -(wr+wl);

if strcmp(bc(1), 'dirichlet')
    wc(1) = 1;
    wr(1) = 0;
end

if strcmp(bc(end), 'dirichlet')
    wc(end) = 1;
    wl(end) = 0;
end

A = spdiags([wr' wc' wl'], [-1 0 1], nx, nx)';

end