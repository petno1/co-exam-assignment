clear
clc


H = [2 0; 0 2];
g = [-2; -5];
A = [1 -1 -1 1 0; -2 -2 2 0 1];
b = [-2; -6; -2; 0; 0];

[x, fval, exitflag] = quadprog(H, g', A', b);

if exitflag == 1
    fprintf('Optimal solution found: x = [%f; %f]\n', x(1), x(2));
    fprintf('Optimal objective value: %f\n', fval);
else
    fprintf('Optimal solution not found.\n');
end


%% Solve QP using primal active-set algorithm
[x,lambda,Wset,it] = qpsolverActiveSet(H,g,A,b,x0)