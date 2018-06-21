function sol = PNPTimeIntegrator(x, lambda, v, tf)
    dx = diff(x);
    xc = x(1:end-1) + 0.5*dx;
    sol.grid = struct('x', x, 'xc', xc, 'dx', dx);
    sol.options = struct('lambda', lambda, 'tf', tf, 'v', v, ...
        'iter_max', 5, 'tol', 1e-6, 'dtmax', 0.5*lambda^2);

    nx = length(x);
    cp  = 0.5*ones(nx-1,1);
    cm  = 0.5*ones(nx-1,1);
    psi = linspace(-v, v, nx)';
    
    pnp0 = [cp; cm; psi];
    mass = [[speye(nx-1) sparse(nx-1, nx-1) sparse(nx-1, nx)]; ...
        [sparse(nx-1, nx-1) speye(nx-1) sparse(nx-1, nx)]; ...
        [sparse(nx, nx-1) sparse(nx, nx-1) sparse(nx, nx)]];
    opt = odeset('Mass', mass, 'AbsTol', 1e-7, 'RelTol', 1e-4);
    pnp = ode15s(@(t, y) PNPSystem(t, y, sol.grid, sol.options), [0 tf], pnp0, opt);
    
    sol.t   = pnp.x;
    sol.cp  = pnp.y(1:nx-1, :);
    sol.cm  = pnp.y(nx:2*nx-2, :);
    sol.psi = pnp.y(2*nx-1:end, :);
end

function dpnp = PNPSystem(~, pnp, grid, options)
    nx = length(grid.x);
    cp  = pnp(1:nx-1);
    cm  = pnp(nx:2*nx-2);
    psi = pnp(2*nx-1:end);
    
    e = -grad(grid.x, psi);
    
    fp  = getFlux(grid, cp, e,  1, 'linear');
    fm  = getFlux(grid, cm, e, -1, 'linear');
    
    dcp  = -diff(fp)./grid.dx;
    dcm  = -diff(fm)./grid.dx;
    rho  = options.lambda^(-2)*cell2node(grid.x, cp - cm);
    lap  = diff(diff(psi)./grid.dx)./diff(grid.xc);
    dpsi = [psi(1) + options.v; lap + rho(2:end-1); psi(end) - options.v];
    
    dpnp = [dcp; dcm; dpsi];
end

function f = getFlux(grid, c, e, z, varargin)
    if nargin < 5
        method = 'linear';
    else
        method = lower(varargin{1});
    end
    
    fd = [0; -diff(c)./diff(grid.xc); 0];
    fe = z*e.*cell2node(grid.x, c, method);
    f = fd + fe;
    f(1) = 0; f(end) = 0;
end