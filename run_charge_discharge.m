% v = [0.1 1 2 5];
v = 0.01;
sol.c = {};
sol.p = {};
sol.d = {};

delta = 0.1;
beta  = 0.001;
tf = 5;

colors = {'b', 'r', 'k', 'm', 'g'};

if ishandle(1)
    close 1
end
h = figure(1); hold on;

for i=1:length(v)
    sol.c{i} = PNP([-v(i),v(i)], delta, beta, tf);
    sol.p{i} = pb_dirichlet(v(i), delta, sol.c{i}.x);
    sol.d{i} = PNP_discharge(sol.p{i}, delta, beta, tf);
    
    qinf = sol.d{i}.q(1);    
    plot(sol.c{i}.t, sol.c{i}.q / qinf, '-', 'color', colors{i}, ...
        'linewidth', 1.5);
    line_fewer_markers(sol.d{i}.t, 1 - sol.d{i}.q / qinf, 10, ...
        '--o', 'spacing', 'curve', 'color', colors{i}, ...
        'markerfacecolor', colors{i}, 'markersize', 6, 'linewidth', 1.5);
    set(gca, 'fontsize', 14);
    xlabel('$tD^+/L^2$', 'fontsize', 16, 'interpreter', 'latex');
    ylabel('Charging Fraction $(\eta)$', 'fontsize', 16, 'interpreter', 'latex');     
    xlim([0, tf+0.01]);
    ylim([0, 1]);
    axis square;   
    drawnow;        
end

h = tightfig(h);
box on;
    
% saveas(h, sprintf('beta_%1.2f_delta_%1.2f.pdf', beta, delta));

hold off
%%
h = figure(1);

idx = 4;
s = sol.c{idx};
pb = sol.p{idx};
tc = 1:10:length(s.t);
% frames(length(tc)) = struct('cdata',[],'colormap',[]);
fc = 1;
for i=tc
    plot(s.y, s.cp(:,i), 'r', 'linewidth', 1.5); hold on
    plot(s.y, s.cm(:,i), 'b', 'linewidth', 1.5);
    plot(s.y, pb.cp, 'r--', 'linewidth', 0.5);
    plot(s.y, pb.cm, 'b--', 'linewidth', 0.5); hold off
        
    set(gca, 'fontsize', 14);
    title(sprintf('$T = %1.2f$', s.t(i)), 'fontsize', 18, 'interpreter', 'latex');    
    xlabel('$x$', 'interpreter', 'latex', 'fontsize', 20);
    ylabel('$c\pm$', 'interpreter', 'latex', 'fontsize', 20);
    legend({'$c^+$', '$c^-$'}, 'interpreter', 'latex', 'fontsize', 16, ...
        'location', 'southeast');
    ylim([0, 1.5]);
    axis square;
    drawnow;
%     frames(fc) = getframe;
    
    saveas(gcf, sprintf('imgs/c%d_%04d.png', idx, fc));
    
    fc = fc + 1;
end

command = sprintf(strcat('/usr/local/bin/ffmpeg -r 15 -i imgs/c%d_%%04d.png', ...
    ' -vcodec libx264 -pix_fmt yuv420p vids/c%d.mp4'), idx, idx);
system(command)

%%
% colors = {'b', 'r', 'k', 'm', 'g'};
% 
% if ishandle(1)
%     close 1
% end
% h = figure(1); hold on;
% 
% for i=1:length(v)
%     qinf = sol.d{i}.q(1);    
%     plot(sol.c{i}.t, sol.c{i}.q / qinf, '-', 'color', colors{i}, ...
%         'linewidth', 1.5);
%     line_fewer_markers(sol.d{i}.t, 1 - sol.d{i}.q / qinf, 10, ...
%         '--o', 'spacing', 'curve', 'color', colors{i}, ...
%         'markerfacecolor', colors{i}, 'markersize', 6, 'linewidth', 1.5);
%     set(gca, 'fontsize', 14);
%     xlabel('$tD^+/L^2$', 'fontsize', 16, 'interpreter', 'latex');
%     ylabel('Charging Fraction $(\eta)$', 'fontsize', 16, 'interpreter', 'latex');     
%     xlim([0, tf+0.01]);
%     ylim([0, 1]);
%     axis square;   
%     drawnow;       
% end
% 
% h = tightfig(h);
% box on;
% saveas(h, 'beta_0.1_delta_0.1.pdf')
% 
% hold off
