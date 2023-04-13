H = [2, 0; 0, 8];
f = [-8; -16];

% equalities
Aeq = [1, 2];
beq = 12;

% inequalities (converted to greater-equal form)
A = [-1, -1; -1, 0; 0, -1];
b = [-10; -1; -6];

% bounds
lb = [1; 1];
ub = [3; 6];

[n,m] = size(A);

%% Convert bl <= A' x <= bu, l <= x <= u to c(x) = Abar' x + bbar >= 0
Abar = [full(A) full(-A) eye(n,n) -eye(n,n)]; 
bbar = [b; zeros(m,1); -lb; ub];

x0=[2;3];

size(Abar)



%% Solve QP using primal active-set algorithm
[x,lambda,Wset,it] = qpsolverActiveSet(H,g,Abar,bbar,x0);
