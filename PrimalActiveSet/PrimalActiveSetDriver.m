%% Step 1 - prep the inputs
n=400;
[H, g, A, b, C, dl, du, l, u] = randomQPGenerator(n,0.90);
% Convert dl <= C' x <= cu, l <= x <= u to c(x) = Cbar' x + dbar >= 0
Cbar = [C -C eye(n,n) -eye(n,n)];
dbar = [-dl; du; -l; u];
[~,m] = size(Cbar);
[~,j] = size(A);


[x0] = linprog(zeros(1,n)',[A,Cbar]',[b;dbar]);

[x] = PrimalActiveSet(H, g, A, b, Cbar, dbar, x0);
[x2] = quadprog(H, g', Cbar',dbar,A',b);


0.5*x'*H*x-g'*x
0.5*x2'*H*x2-g'*x2