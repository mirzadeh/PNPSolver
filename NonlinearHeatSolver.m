% solve dc/dt = div(exp(c) grad c) with c(1) = 1 and c(end) = 0 as boundary
% conditions
x = linspace(0, 1)';
t = linspace(0, 0.2);

s = nonlinearSolver(x, t);

figure(1); clf; hold on
plot(s.t, s.c(1:15:end, :), 'k.-', 'LineWidth', 1.5, ...
    'MarkerSize', 12, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'r');

function sol = nonlinearSolver(x, t)
    dx = x(2) - x(1);
    c0 = x.*(1-x);
    opt = odeset('AbsTol', 1e-6, 'RelTol', 1e-3);
    ode = ode15s(@odesystem, [0 t(end)], c0, opt);
    sol.t = ode.x;
    sol.c = ode.y;
    
    function dc = odesystem(~, c)
        k = exp(2*(c(1:end-1) + 0.5*diff(c)).^2) + 0.01;
        f = [0; k.*diff(c); 0]/dx;
        dc = diff(f)/dx;
    end
end