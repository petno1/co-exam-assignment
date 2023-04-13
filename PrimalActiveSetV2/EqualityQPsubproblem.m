function [x, lagrangeMultipliers] = EqualityQPsubproblem(H,g,A,b,Aeq,beq)
% Combine inequality and equality constraints
A = [A Aeq];
b = [b; beq];
[n,m] = size(A);
% Construct the KKT system
K = [H, A';
     A, zeros(m,m)];
d = [-g; b];
% Solve the KKT system
size(K)
size(d)
z = K \ d;
% Extract the solution
x = z(1:size(H,1));
lagrangeMultipliers = z(size(H,1)+1:end);
end
