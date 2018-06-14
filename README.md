![](figures/cover.png)
This is a simple PNPSolver based on hybrid FD/FV formulation.

* PNPSolver: Solves PNP equation between two blocking electrodes.
Useage:
``` MATLAB
>> x = gridGen(-1, 1, 100, 1.5); % stretch grid near boundaries
>> pnp = PNPSolver(x, 0.1, 3, 2); % solve pnp with lambda = 0.1 (edl/L), V = 3 (KT/e), and tf = 2 (L^2/D)
>> plot(pnp.grid.x, pnp.psi(:,end), pnp.grid.xc, pnp.cp(:,end) - pnp.cm(:,end), 'r'); % plot final solution
```

* PBSolver: Solves PB equation between two electrodes of fixed potential
Useage:
``` MATLAB
>> x = gridGen(-1, 1, 100, 1.5);
>> pb = PBSolver(x, 0.1, 3);
>> plot(pb.grid.x, pb.psi);
```
* MembraneSolver: Solves PNP equation between two electrodes and a perm-selective membrane in
  between
