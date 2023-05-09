function [x_sol] = InteriorPointQP(H, g, A, b, C, d, x0)

% Convert Cx >= d into an equality constraint using a slack variable
%   A' x = b
%   C' x + d + s = 0
%   s > 0

% Sets constants for the algorithm
mIn = length(x0);
tol = 0.000001;
eta = 0.99;

y0 = zeros(size(b));
z0 = ones(size(d));
s0 = ones(size(d));

x=x0;
y=y0; %(lambda)
z=z0; %(mu)
s=s0;

S = diag(s);
Z = diag(z);
e = ones(size(s));

rL = H*x+g-A*y-C*z;
rA = -A'*x+b;
rC = -C'*x+s+d;
rSZ = S*Z*e;

dualGap = (z'*s)/(mIn);

% Converged
Converged = (norm(rL,inf) <= tol) && ...
            (norm(rA,inf) <= tol) && ...
            (norm(rC,inf) <= tol) && ...
            (abs(s) <= tol);
%%       
max_iter = 100;
iter = 0;
hist = [];

while ~Converged && (iter<max_iter)
    iter = iter+1;
 
    disp('===================')
    iter
    disp('===================')

    h_bar = H + C*(Z/S)*C';
    rl_bar =  rL-C*(S\Z)*(rC-Z\rSZ);
    KKT = [h_bar -A; -A' zeros(size(A,2))];
    [L,D,p] = ldl(KKT,"lower",'vector');    
    rhs = -[rl_bar;rA];
    solution(p) = L'\(D\(L\rhs(p)));

    %% Affine Direction
    dx_affine = solution(1:length(x))';
    dz_affine =-(S\Z)*C'*dx_affine+(S\Z)*(rC-Z\rSZ);
    ds_affine = -(Z\rSZ)-(Z\(S*dz_affine));

    dZS = [dz_affine; ds_affine];
    alphas = (-[z; s] ./ dZS);
    alphaMax_affine = eta*min([1; alphas(dZS < 0)]);

    %% Duality gap and centering parameter
    affineDualGap = ((z+alphaMax_affine*dz_affine)'*(s+alphaMax_affine*ds_affine))/(mIn);
    sigma = (affineDualGap./dualGap)^3;

    %% Affine-Centering-Correction Direction
    rSZ_bar = rSZ + (ds_affine.*dz_affine.*e) - (affineDualGap*sigma.*e);
    rl_bar = rL- C*(S\Z)*(rC-Z\rSZ_bar);

    [L,D,p] = ldl(KKT,'vector'); 

    rhs = -[rl_bar ; rA];
    solution(p) = L'\(D\(L\rhs(p)));

    dx = solution(1:length(x))';
    dy = solution(length(x)+1:length(x)+length(y))';
    
    dz = -(S\Z)*C'*dx+(S\Z)*(rC-Z\rSZ_bar);
    ds = -(Z\rSZ_bar)-(Z\(S*dz));

    dZS = [dz;ds];
    alphas =(-[z;s]./dZS);
    alpha = min([1; alphas(dZS<0)])*eta;

    % Update of position
    x = x + alpha * dx;
    y = y + alpha * dy;
    z = z + alpha * dz;
    s = s + alpha * ds;

    S = diag(s);
    Z = diag(z);
    e = ones(size(s));

    rL = H*x+g-A*y-C*z;
    rA = -A'*x+b;
    rC = -C'*x+s+d;
    rSZ = S*Z*e;

    dualGap = (z'*s)/(mIn);

    Converged = (norm(rL,inf) <= tol) && ...
            (norm(rA,inf) <= tol) && ...
            (norm(rC,inf) <= tol) && ...
            (abs(s) <= tol);

hist = vertcat(hist, 0.5*x'*H*x-g'*x);

end

hist
iter
x_sol = x;

end

