%% Step 1 - prep the inputs
n=30;
[H, g, A, b, C, dl, du, l, u] = randomQPGenerator(n,0.20);
% Convert dl <= C' x <= cu, l <= x <= u to c(x) = Cbar' x + dbar >= 0
Cbar = [C -C eye(n,n) -eye(n,n)];
dbar = [-dl; du; -l; u];
[~,m] = size(Cbar);
[~,j] = size(A);

[x] = InteriorPointQP(H, g, A, b, Cbar, dbar, zeros(n,1));
[x2] = quadprog(H, g', Cbar',dbar,A',b);

0.5*x'*H*x-g'*x
0.5*x2'*H*x2-g'*x2
