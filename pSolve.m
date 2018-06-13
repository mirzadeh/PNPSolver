function pn = pSolve(x,cp,cm,kappa,bc)
% interpolate charge density on nodes
rho = kappa^2*(cp-cm);
A = matGen(x);
f = -cell2node(x, rho);

% adjust boundary conditions
f(1) = bc(1);
f(end) = bc(end);

pn = A\f;
end