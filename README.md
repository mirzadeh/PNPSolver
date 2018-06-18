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
>> x.l = gridGen(-1, 0, 50, 1.5);
>> x.r = gridGen( 0, 1, 50, 1.5);
>> v = 1; lambda = 0.1; lambda_m = 0.1; G = 10; tf = 1;
>> mem = MembraneSolve(x, v, G, lambda, lambda_m, tf);
>> plot(mem.t, mem.cp.r(1,:) - mem.cp.l(end,:)); % plot \delta c^+ across membrane
```
![](figures/dcp.png)
For more help take a look at [runMembrane](runMembrane.m) script.
