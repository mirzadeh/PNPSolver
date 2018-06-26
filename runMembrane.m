% grid
x.l = gridGen(-1, 0, 50, 1.5);
x.r = gridGen( 0, 1, 50, 1.5);

% physical parameters
options = struct('v', 1, 'lambda', 0.1, 'lambda_m', 0.1, 'tf', 1);
ions = {ion(1, 0.5, 1) ion(-1, 0.5)};

sol = MembraneSolver(x, options, ions);
%% plot
figure(1); clf; hold on;
plot(sol.t, sol.mem.dc(1, :), 'r', sol.t, sol.mem.dc(2, :), 'b', ...
     sol.t, sol.mem.dpsi, 'k', 'LineWidth', 1.5);

xlabel('$t D/L^2$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$\Delta c^\pm, \Delta \psi$', 'Interpreter', 'latex', 'FontSize', 12);
legend({'$\Delta c^+$', '$\Delta c^-$', '$\Delta \psi$'}, ...
    'Interpreter', 'latex', 'FontSize', 12);
set(gca, 'FontSize', 12);

print -f1 -dpng -r200 figures/delta.png

figure(2); clf; hold on
plot(sol.t, sol.mem.j(1, :), 'r', sol.t, sol.mem.j(2, :), 'b', 'LineWidth', 1.5);

xlabel('$t D/L^2$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$j_{mem}^\pm$', 'Interpreter', 'latex', 'FontSize', 12);
legend({'$j_{mem}^+$', '$j_{mem}^-$'}, 'Interpreter', 'latex', ...
    'FontSize', 12, 'Location', 'southeast');
set(gca, 'FontSize', 12);

print -f2 -dpng -r200 figures/jmem.png

figure(3); clf; hold on
plot(sol.grid.xc.l, sol.ions{1}.c.l(:,end), 'r', 'LineWidth', 1.5);
plot(sol.grid.xc.l, sol.ions{2}.c.l(:,end), 'b', 'LineWidth', 1.5);
plot(sol.grid.x.l, sol.psi.l(:,end), 'k', 'LineWidth', 1.5);

plot(sol.grid.xc.r, sol.ions{1}.c.r(:,end), 'r', 'LineWidth', 1.5);
plot(sol.grid.xc.r, sol.ions{2}.c.r(:,end), 'b', 'LineWidth', 1.5);
plot(sol.grid.x.r, sol.psi.r(:,end), 'k', 'LineWidth', 1.5);

xlabel('$x/L$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$c^\pm, \psi$', 'Interpreter', 'latex', 'FontSize', 12);
legend({'$c^+$', '$c^-$', '$\psi$'}, 'Interpreter', 'latex', ...
    'FontSize', 12, 'Location', 'southeast');
set(gca, 'FontSize', 12);

print -f3 -dpng -r200 figures/spatial.png


