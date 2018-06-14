function pn = pSolve(x, cp, cm, kappa, bc)
% interpolate charge density on nodes
f = cell2node(x, kappa^2*(cp-cm));

% adjust boundary conditions
f(1) = bc(1);
f(end) = bc(end);

pn = matGen(x) \ f;
end