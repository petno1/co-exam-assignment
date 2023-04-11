function [x, lagrangeMultipliers] = EqualityQPsubproblem(H,g,A,b)

[n,m] = size(A);
K = [H -A; -A' zeros(m,m)];
d = [-g; b];
z = K\d;
x = z(1:n,1);
lagrangeMultipliers = z(n+1:n+m,1);
