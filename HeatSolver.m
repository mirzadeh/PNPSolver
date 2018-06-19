x = linspace(0, 1)';
t = linspace(0, 1);

s1 = backwardEuler(x, t);
s2 = matlabIntegrator(x, t);

figure(1); clf; hold on
plot(s1.t, s1.c(20:20:80, :), 'k', 'LineWidth', 1.5);
plot(s2.t, s2.c(20:20:80, :), 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'w');

function sol = backwardEuler(x, t)
    A = matGen(x);
    dt = t(2) - t(1);
    D = 1/dt*speye(length(x));
    D(1, 1) = 0; D(end, end) = 0;
    
    sol.c = zeros(length(x), length(t));
    sol.c(1,:) = 1;
    sol.t = 0;
    for i = 1:length(t)-1
        f = sol.c(:, i)/dt;
        f(1) = 1; f(end) = 0;
        sol.c(:, i+1) = (A + D) \ f;
        sol.t = cat(1, sol.t, t(i+1));
    end
end

function sol = matlabIntegrator(x, t)
    A = matGen(x);
    c0 = zeros(size(x));
    c0(1) = 1;
    opt = odeset('AbsTol', 1e-6, 'RelTol', 1e-2);
    ode = ode15s(@(t, c) odesystem(t, c, A), [0 t(end)], c0, opt);
    sol.t = ode.x;
    sol.c = ode.y;
    
    function dc = odesystem(~, c, A)
        dc = -A*c;
        
        % apply bc
        dc(1) = 0;
        dc(end) = 0;
    end
end