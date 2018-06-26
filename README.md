![](figures/cover.png)

This is a simple PNPSolver based on hybrid FD/FV formulation.

* **PNPSolver**: Solves PNP equation between two blocking electrodes. Useage:
``` MATLAB
>> x = gridGen(-1, 1, 100, 1.5); 
>> lambda = 0.1; v = 3; tf = 2;
>> pnp = PNPSolver(x, lambda, v, tf);
>> plot(pnp.grid.x, pnp.psi(:,end), pnp.grid.xc, pnp.cp(:,end) - pnp.cm(:,end), 'r'); 
>> mesh(pnp.t, pnp.grid.xc, pnp.cp);
```

* **PBSolver**: Solves PB equation between two electrodes of fixed potential. Useage:
``` MATLAB
>> x = gridGen(-1, 1, 100, 1.5);
>> lambda = 0.1; v = 3;
>> pb = PBSolver(x, lambda, v);
>> plot(pb.grid.x, pb.psi);
```
* **MembraneSolver**: Solves PNP equation between two electrodes and a perm-selective membrane in between. Usage:
``` MATLAB
>> x.l = gridGen(-1, 0, 50, 1.5); % left grid
>> x.r = gridGen( 0, 1, 50, 1.5); % right grid
>> options = struct('v', 1, 'lambda', 0.1, 'lambda_m', 0.1, 'tf', 1); % solver options
>> ions = {ion(1, 0.5, 1) ion(-1, 0.5)}; % a binary electrolyte with G = 1 for cations
>> sol = MembraneSolve(x, options, ions);
>> plot(sol.t, sol.mem.j(1,:), 'r', sol.t, sol.mem.j(2,:), 'b'); % plot j_{mem}^\pm across membrane
```
![](figures/jmem.png)
For more help take a look at [runMembrane](runMembrane.m) script.
