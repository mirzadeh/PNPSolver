function sol = MembraneSolver(x, options, ions)
    sol.options = options;
    
    % default to binary symmetric electrolyte if no ions are supplied
    if nargin < 3
        ions{1} = struct('z',  1, 'c0', 0.5, 'd', 1, 'G', 1);
        ions{2} = struct('z', -1, 'c0', 0.5, 'd', 1, 'G', 0);
    end
    sol.ions = ions;

    % prepare the grid
    dx = struct('l', diff(x.l), 'r', diff(x.r));
    xc.l = x.l(1:end-1) + 0.5*dx.l; 
    xc.r = x.r(1:end-1) + 0.5*dx.r;    
    sol.grid = struct('x', x, 'xc', xc, 'dx', dx);
        
    % compute the initial conditions        
    nl = length(x.l); nr = length(x.r);
    c0 = zeros(nl+nr-2, length(ions));
    p0 = zeros(nl+nr, 1);
    for i = 1:length(ions)
        c0(:, i) = ions{i}.c0*ones(nl+nr-2, 1);
    end
    p0(1:nl) = linspace(-options.v, 0, nl);
    p0(nl+1:end) = linspace(0, options.v, nr);

    % integrate the PNP equations -- note that the PNP equations is a index-1
    % DAE since the Poisson eqution does not involve time derivatives.
    % Therefore we need to define the mass matrix
    nc = nl+nr-2; np = nl+nr;
    dofs = nc*length(ions) + np;
    diag = ones(dofs, 1);
    diag(end-np+1:end) = 0;        
    mass = spdiags(diag, 0, dofs, dofs);
    
    % set the options for ODE solver and call it
    opt = odeset('Mass', mass, 'OutputFcn', @PNPPrint, 'AbsTol', 1e-6, 'RelTol', 1e-3);
    pnp = ode15s(@(t, y) PNPSystem(t, y, sol.grid, options, ions), ...
        [0 options.tf], [reshape(c0, [], 1); p0], opt);
    
    % extract solutions
    sol.t = pnp.x';
    
    sol.psi.l = pnp.y(end-np+1:end-nr, :);
    sol.psi.r = pnp.y(end-nr+1:end, :);
    
    % membrane variables
    sol.mem.dpsi = squeeze(sol.psi.r(1,:) - sol.psi.l(end,:));
    sol.mem.dc = zeros(length(ions), length(sol.t));
    sol.mem.j = zeros(length(ions), length(sol.t));
    
    cn = reshape(pnp.y(1:nc*length(ions), :), nc, length(ions), []);    
    for i = 1:length(ions)
        sol.ions{i}.c.l = squeeze(cn(1:nl-1, i, :));
        sol.ions{i}.c.r = squeeze(cn(nl:nl+nr-2, i, :));
        
        % store membrane values for easy use        
        cl = squeeze(interp1(sol.grid.xc.l, sol.ions{i}.c.l, 0, 'linear', 'extrap'));
        cr = squeeze(interp1(sol.grid.xc.r, sol.ions{i}.c.r, 0, 'linear', 'extrap'));
        sol.mem.dc(i,:) = cr - cl;
        sol.mem.j(i,:) = -ions{i}.G*(sol.mem.dpsi + log(cr./cl));
    end        
end

function r = PNPSystem(~, y, grid, options, ions)
    % Compute the residual of the PNP equations. This is the fully implicit
    % formulation. We store everything in a single vector which includes all
    % ionic concentrations followed by the potential            
    nl = length(grid.x.l);
    nr = length(grid.x.r);
    nc = nl+nr-2;
    np = nl+nr;
    
    dc = zeros(nc, length(ions));
    z = zeros(1, length(ions));
    
    psi = y(end-np+1:end);
    c = reshape(y(1:nc*length(ions)), nc, length(ions));
    
    el = -grad(grid.x.l, psi(1:nl));
    er = -grad(grid.x.r, psi(nl+1:end));
    
    for i = 1:length(ions)
        z(i) = ions{i}.z;
        fl = getFlux(grid.x.l, c(1:nl-1, i), el, ions{i}, 'linear');
        fr = getFlux(grid.x.r, c(nl:end, i), er, ions{i}, 'linear');
        
        % compute membrane flux
        cl = interp1(grid.xc.l, c(1:nl-1, i), 0, 'linear', 'extrap');
        cr = interp1(grid.xc.r, c(nl:end, i), 0, 'linear', 'extrap');
        jmem = -ions{i}.G*(psi(nl+1) - psi(nl) + log(cr/cl));
        
        % apply boundary conditions
        fl(1) = 0; fl(end) = jmem;
        fr(1) = jmem; fr(end) = 0;
        
        dc(:, i) = [-diff(fl)./grid.dx.l; -diff(fr)./grid.dx.r];
    end
    
    % Poisson equation -- Note we impose constant potential boundary conditions
    rhol = options.lambda^(-2)*cell2node(grid.x.l, sum(z.*c(1:nl-1,:), 2));
    rhor = options.lambda^(-2)*cell2node(grid.x.r, sum(z.*c(nl:end,:), 2));    
    dpsi = [laplace(grid.x.l, psi(1:nl)) + rhol; laplace(grid.x.r, psi(nl+1:end)) + rhor];
    
    % apply electrode boundary conditions
    dpsi([1, end]) = [psi(1) + options.v; psi(end) - options.v];
    
    % apply membrane boundary conditions
    dpsi(nl) = psi(nl+1) - psi(nl) - options.lambda_m*(psi(nl)-psi(nl-1))/grid.dx.l(end);
    dpsi(nl+1) = (psi(nl+2) - psi(nl+1))/grid.dx.r(1) - (psi(nl) - psi(nl-1))/grid.dx.l(end);
    
    % pack everything
    r = [reshape(dc, [], 1); dpsi];
    
    function f = getFlux(x, c, e, ion, varargin)
        if nargin < 5
            method = 'linear';
        else
            method = lower(varargin{1});
        end
        
        fd = -grad(x, c, 'cell');
        fe = ion.z*e.*cell2node(x, c, method);
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
