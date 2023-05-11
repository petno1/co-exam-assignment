 
%% Step 1 - prep the inputs
n=500;
[H, g, A, b, C, dl, du, l, u] = randomQPGenerator(n,0.90);
%% Convert dl <= C' x <= cu, l <= x <= u to c(x) = Cbar' x + dbar >= 0
Cbar = [C -C eye(n,n) -eye(n,n)];
dbar = [-dl; du; -l; u];


opt = mpcInteriorPointOptions;

[x1] = InteriorPointQP(H, g, A, b, Cbar, dbar, zeros(n,1));
[x2] = mpcInteriorPointSolver(H, g, full(Cbar)',dbar,A',b,zeros(n,1),opt);

0.5*x1'*H*x1-g'*x1
0.5*x2'*H*x2-g'*x2