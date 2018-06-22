function sol = PNPTimeIntegrator(x, lambda, v, tf)
    dx = diff(x);
    xc = x(1:end-1) + 0.5*dx;
    sol.grid = struct('x', x, 'xc', xc, 'dx', dx);
    sol.options = struct('lambda', lambda, 'v', v, 'tf', tf);
    
    % compute the initial conditions
    nx = length(x);
    cp  = 0.5*ones(nx-1,1);
    cm  = 0.5*ones(nx-1,1);
    psi = linspace(-v, v, nx)';
    
    % integrate the PNP equations -- note that the PNP equations is a index-1
    % DAE since the Poisson eqution does not involve time derivatives.
    % Therefore we need to define the mass matrix
    mass = [[speye(nx-1) sparse(nx-1, nx-1) sparse(nx-1, nx)]; ...
        [sparse(nx-1, nx-1) speye(nx-1) sparse(nx-1, nx)]; ...
        [sparse(nx, nx-1) sparse(nx, nx-1) sparse(nx, nx)]];
    opt = odeset('Mass', mass, 'OutputFcn', @PNPPrint);
    pnp = ode15s(@(t, y) PNPSystem(t, y, sol.grid, sol.options), [0 tf], [cp; cm; psi], opt);
    
    sol.t   = pnp.x;
    sol.cp  = pnp.y(1:nx-1, :);
    sol.cm  = pnp.y(nx:2*nx-2, :);
    sol.psi = pnp.y(2*nx-1:end, :);
end

function dpnp = PNPSystem(~, pnp, grid, options)
    % Compute the residual of the PNP equations. This is the fully implicit
    % formulation. We store everything in a single vector which includes cp,
    % cm, and psi.
    nx = length(grid.x);
    cp  = pnp(1:nx-1);
    cm  = pnp(nx:2*nx-2);
    psi = pnp(2*nx-1:end);    
    e = -grad(grid.x, psi);
    
    % Nernst-Planck equations -- Note we impose zero flux boundary conditions
    % on both ends
    % cp
    fp  = getFlux(grid, cp, e,  1, 'linear'); 
    fp(1) = 0; fp(end) = 0;
    dcp = -diff(fp)./grid.dx;
    
    % cm
    fm  = getFlux(grid, cm, e, -1, 'linear');    
    fm(1) = 0; fm(end) = 0;    
    dcm = -diff(fm)./grid.dx;
    
    % Poisson equation -- Note we impose constant potential boundary conditions
    rho  = options.lambda^(-2)*cell2node(grid.x, cp - cm);
    dpsi = laplace(grid.x, psi) + rho;   
    dpsi([1, end]) = [psi(1) + options.v; psi(end) - options.v];
    
    % pack everything
    dpnp = [dcp; dcm; dpsi];
    
    function f = getFlux(grid, c, e, z, varargin)
        if nargin < 5
            method = 'linear';
        else
            method = lower(varargin{1});
        end
        
        fd = -grad(grid.x, c, 'cell');
        fe = z*e.*cell2node(grid.x, c, method);
        f = fd + fe;        
    end
end

function status = PNPPrint(t, ~, flag)    
    switch flag
        case 'done'
            fprintf('\n');
        otherwise
            fprintf('t = %f\n', t(1));
    end
    status = 0;
end
