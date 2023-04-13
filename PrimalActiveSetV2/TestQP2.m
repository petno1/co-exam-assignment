clear
clc


%% Parameters for random QP
n=20; 
density=0.5; 

%% Initial guess of solution (we know zero is feasible)
x0 = zeros(n,1);

%% Generate random QP
[H, g, A, b, C, dl, du, l, u] =randomQPGenerator(n, density); 

%% Convert bl <= A' x <= bu, l <= x <= u to c(x) = Abar' x + bbar >= 0
Cbar = [full(C) full(-C) eye(n,n) -eye(n,n)]; 
dbar = [-dl; du; -l; u];

%% Solve QP using primal active-set algorithm

[x1] = qpsolverActiveSet(H,g,Cbar,dbar,A,b,x0)
[x2] = quadprog(H,g',Cbar',dbar,A',b)