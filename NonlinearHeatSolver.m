% solve dc/dt = div(exp(c) grad c) with c(1) = 1 and c(end) = 0 as boundary
% conditions
x = linspace(0, 1)';
t = linspace(0, 1);

s = nonlinearSolver(x, t);

figure(1); clf; hold on
plot(s.t, s.c(20:20:80, :), 'k.-', 'LineWidth', 1.5, ...
    'MarkerSize', 12, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'r');
ylim([0, 1]);

function sol = nonlinearSolver(x, t)
    dx = x(2) - x(1);
    c0 = zeros(size(x));
    c0(1) = 1;
    opt = odeset('AbsTol', 1e-6, 'RelTol', 1e-2);
    ode = ode15s(@odesystem, [0 t(end)], c0, opt);
    sol.t = ode.x;
    sol.c = ode.y;
    
    function dc = odesystem(~, c)
        dc = zeros(size(c));
        k = exp(c(1:end-1) + 0.5*diff(c).^2) + 0.01;
        dc(2:end-1) = diff(k.*diff(c))/dx^2;
    end
end