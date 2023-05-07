function [x_sol] = InteriorPointQP(H, g, A, b, C, d, x0)

% Convert Cx >= d into an equality constraint using a slack variable
%   A' x = b
%   C' x + d + s = 0
%   s > 0


% Sets constants for the algorithm
mIn = length(d);
m = length(b);
n = length(x0);
tol = 0.000001;
eta = 0.995;

y0 = zeros(size(b));
z0 = ones(size(d));
s0 = ones(size(d));


x=x0;
y=y0; %(lambda)
z=z0; %(mu)
s=s0;


rL = H*x+g-A*y-C*z;
rA = -A'*x+b;
rC = -C'*x+s+d;

S = diag(s);
M = diag(z);
e = ones(size(s));

rSm = S*M*e;

% Converged
Converged = (norm(rL,inf) <= tol) && ...
            (norm(rA,inf) <= tol) && ...
            (norm(rC,inf) <= tol) && ...
            (abs(s) <= tol);
%%       
max_iter = 100;
iter = 0;

dx = -rL;
dy = -rA;
dz = -rC;
ds = -rSm;

while ~Converged && (iter<max_iter)
    iter = iter+1;
 
    rl_bar = rL - C*(M/S)*(rC-inv(M)*rSm);
    h_bar = H + C*(M/S)*C';
    KKT = [h_bar -A; -A' zeros(m)];
    [L,D,~] = ldl(KKT,'vector');    
    rhs = -[rl_bar ; rA ];
    solution = L'\(D\(L\rhs));
    dx = solution ( 1 : length (x) )';
    dy = solution (length(x)+1:length(x)+length(y))'; 
    dz = (M/S)*C'.*dx + (M/S)*(rC-inv(M)*rSm);
    ds =  (inv(M)*rSm)-(inv(M)*S*dz);
    
    dZS = [dz; ds];
    alphas = (-[z; s] ./ dZS);
    alpha = min([1; alphas(dZS < 0)]);
    alphaBar = eta * alpha;
    
    % Update of position
    x = x + alphaBar * dx
    y = y + alphaBar * dy
    z = z + alphaBar * dz
    s = s + alphaBar * ds


    




end




x_sol = 0

end

