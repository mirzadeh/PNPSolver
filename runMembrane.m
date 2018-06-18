% grid
x.l = gridGen(-1, 0, 50, 2);
x.r = gridGen( 0, 1, 50, 2);

% physical parameters
v = 1;
lambda = 0.1; % edl thickness
lambda_m = 0.1;
G = [0 0.1 1 10];
tf = 1;

sol = cell(length(G), 1);
for i = 1:length(G)
    sol{i} = MembraneSolver(x, v, G(i), lambda, lambda_m, tf);
end

%% plot
close all;
for i = 1:length(G)
    figure(1); hold on;
    tplot(sol{i}.t, sol{i}.cp.r(1,:) - sol{i}.cp.l(end,:), 'r', '$\Delta c^+$');
    
    figure(2); hold on;
    tplot(sol{i}.t, sol{i}.cm.r(1,:) - sol{i}.cm.l(end,:), 'b', '$\Delta c^-$');    
    
    figure(3); hold on;
    tplot(sol{i}.t, sol{i}.psi.r(1,:) - sol{i}.psi.l(end,:), 'k', '$\Delta \phi$');    
end

print -f1 -dpng -r100 'figures/dcp.png'
print -f2 -dpng -r100 'figures/dcm.png'
print -f3 -dpng -r100 'figures/dphi.png'

function tplot(x, y, c, ylabel_)
plot(x, y, c, 'LineWidth', 2);
xlabel('$t D/L^2$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel(ylabel_, 'Interpreter', 'latex', 'FontSize', 18);
set(gca, 'FontSize', 16);
end