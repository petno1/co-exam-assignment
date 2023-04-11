% 02612 Constrained Optimization
close all
clear
clc

%% step 1
% test the random qp generator using quadprog
% min φ(x) = 1/2 x′Hx + g′x 
% s.t. A′x + b = 0 
% dl ≤ C′x ≤ du 
% l ≤ x ≤ u 
n=100;
[H, g, A, b, C, dl, du, l, u] = randomQPGenerator(n,0.5);
%% Convert dl <= C' x <= cu, l <= x <= u to c(x) = Cbar' x + dbar >= 0
Cbar = [C -C eye(n,n) -eye(n,n)];
dbar = [-dl; du; -l; u];

%% step 2
% find a well defined problem in the form
% min φ(x) = 1/2 x′Hx + g′x 
% s.t. A′x + b = 0 
% dl ≤ C′x ≤ du 
% l ≤ x ≤ u 

[x] = quadprog(H, g', Cbar', dbar, A', b);


