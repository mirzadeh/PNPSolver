b =  [0 1e-3 1e-2 1e-1 0.5 1];
sol = {};

for i=1:length(b)
    sol{i} = PNP([-1,1], 0.01, b(i), 10);
end

%% plot
qmax = max(sol{end}.q);
close 1;
figure(1); hold on;

for i=1:length(sol)
    plot(sol{i}.t, sol{i}.q/qmax, 'linewidth', 1.5); 
end
axis square

xlabel('$tD^+/L^2$', 'fontsize', 20, 'interpreter', 'latex');
ylabel('$q/q_\infty$', 'fontsize', 20, 'interpreter', 'latex');
ylim([0,1]);
xlim([0,sol{1}.runData.tf]);
set(gca, 'fontsize', 12);

savefig('charge_discharge','pdf');

hold off;