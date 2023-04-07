clear
clc
%% Parameters for random QP
n=10; 
m=20;
alpha=0.1; 
density=0.75; 
%% Initial guess of solution (we know zero is feasible)
x0 = zeros(n,1);


%% Generate random QP
[H, g, A, b, C, dl,du, l, u] = RandomQPGenerator(n,m,alpha,density);

d = dl-du;

[x] = quadprog(H,g',C',d,A',b,l,u,[],optimset('Display','iter'));
