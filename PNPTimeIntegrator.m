function sol = PNPTimeIntegrator(x, options, varargin)
    sol.options = options;
    
    % default to binary symmetric electrolyte if no ions are supplied
    if nargin < 3
        ions{1} = struct('z',  1, 'c0', 0.5, 'd', 1);
        ions{2} = struct('z', -1, 'c0', 0.5, 'd', 1);
    else
        ions = varargin{1};
    end
    sol.ions = ions;

    % prepare the grid
    dx = diff(x);
    xc = x(1:end-1) + 0.5*dx;
    sol.grid = struct('x', x, 'xc', xc, 'dx', dx);
        
    % compute the initial conditions
    nx = length(x);
    y0 = zeros(length(ions)*(nx-1) + nx, 1);
    
    for i = 1:length(ions)
        y0(1 + (i-1)*(nx-1) : i*(nx-1)) = ions{i}.c0*ones(nx-1, 1);
    end
    y0(end-nx+1:end) = options.v*linspace(-1, 1, nx)';
    
    % integrate the PNP equations -- note that the PNP equations is a index-1
    % DAE since the Poisson eqution does not involve time derivatives.
    % Therefore we need to define the mass matrix
    dofs = length(ions)*(nx-1) + nx;
    diag = ones(dofs, 1);
    diag(end-nx+1:end) = 0;
    mass = spdiags(diag, 0, dofs, dofs);
    opt = odeset('Mass', mass, 'OutputFcn', @PNPPrint);
    pnp = ode15s(@(t, y) PNPSystem(t, y, sol.grid, options, ions), [0 options.tf], y0, opt);
    
    % extract solutions
    sol.grid.t = pnp.x';
    for i = 1:length(ions)
        sol.ions{i}.c = pnp.y(1 + (i-1)*(nx-1):i*(nx-1), :);
    end
    sol.psi = pnp.y(end-nx+1:end, :);
end

function dpnp = PNPSystem(~, pnp, grid, options, ions)
    % Compute the residual of the PNP equations. This is the fully implicit
    % formulation. We store everything in a single vector which includes all
    % ionic concentrations followed by the potential
    nx = length(grid.x);
    dc = zeros(nx-1, length(ions));
    zs = zeros(1, length(ions));
    
    psi = pnp(end-nx+1:end);    
    c = reshape(pnp(1:end-nx), nx-1, length(ions));
    
    e = -grad(grid.x, psi);
    
    for i = 1:length(ions)
        zs(i) = ions{i}.z;
        f = getFlux(grid, c(:,i), e, ions{i}, 'linear');
        f(1) = 0; f(end) = 0;
        
        dc(:, i) = -diff(f) ./ grid.dx;
    end    
        
    % Poisson equation -- Note we impose constant potential boundary conditions
    rho  = options.lambda^(-2)*cell2node(grid.x, sum(zs.*c, 2));
    dpsi = laplace(grid.x, psi) + rho;
    dpsi([1, end]) = [psi(1) + options.v; psi(end) - options.v];
    
    % pack everything
    dpnp = [reshape(dc, [], 1); dpsi];
    
    function f = getFlux(grid, c, e, ion, varargin)
        if nargin < 5
            method = 'linear';
        else
            method = lower(varargin{1});
        end
        
        fd = -grad(grid.x, c, 'cell');
        fe = ion.z*e.*cell2node(grid.x, c, method);
        f = ion.d*(fd + fe);
    end
end

function status = PNPPrint(t, ~, flag)
    % a print function to be called after each step of ode solver 
    switch flag
        case 'done'
            fprintf('\n');
        otherwise
            fprintf('t = %f\n', t(1));
    end
    status = 0;
end
