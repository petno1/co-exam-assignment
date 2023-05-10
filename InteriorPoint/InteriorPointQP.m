function [x_sol] = InteriorPointQP(H, g, A, b, C, d, x0)

% Convert Cx >= d into an equality constraint using a slack variable
%   A' x = b
%   C' x + d + s = 0
%   s > 0

% Sets constants for the algorithm
m = length(d);
tol = 1.0e-4;
eta = 0.95;

y0 = zeros(size(b));
z0 = ones(size(d));
s0 = ones(size(d));

S = diag(s0);
Z = diag(z0);
e = ones(size(s0));

x=x0;
y=y0;
z=z0;
s=s0;

rL = H*x+g-A*y-C*z;
rA = -A'*x+b;
rC = -C'*x+s+d;
rSZ = S*Z*e;
dualGap = (z'*s)./m;

% Converged
Converged = (norm(rL) <= tol) && ...
            (norm(rA) <= tol) && ...
            (norm(rC) <= tol) && ...
            (abs(dualGap) <= tol) && ...
            (norm(s) <= tol);
%%       
max_iter = 100;
iter = 0;
hist = [];

hist = vertcat(0.5*x0'*H*x0-g'*x0);


while ~Converged && (iter<max_iter)

    iter = iter+1;

    h_bar = H + C*(Z/S)*C';
    rl_bar =  rL-C*(S\Z)*(rC-Z\rSZ);
    disp('===========')
    KKT = [h_bar -A; -A' zeros(size(A,2))];
    [L,D,p] = ldl(KKT,"lower",'vector');
    rhs = -[rl_bar;rA];

    solution(p) = L'\(D\(L\rhs(p)));

    %% Affine Direction
    dx_affine = solution(1:length(x))';
    dz_affine =-(S\Z)*C'*dx_affine+(S\Z)*(rC-Z\rSZ);
    ds_affine = -(Z\rSZ)-(Z\(S*dz_affine));

    %% Calculate alpha
    dZS = [ dz_affine ; ds_affine ] ;
    alphas = -([z; s]./dZS);
    affineAlpha = min([1;alphas(dZS<0)]);

    %% Duality gap and centering parameter
    affineDualGap = (((z+affineAlpha*dz_affine)'*(s+affineAlpha*ds_affine)))./(m);
    sigma = (affineDualGap/dualGap)^3;

    %% Affine-Centering-Correction Direction
    rSZ_bar = rSZ + (diag(ds_affine)*diag(dz_affine)*e) - (sigma*dualGap*e);
    rl_bar = rL- C*(S\Z)*(rC-Z\rSZ_bar);


    rhs2 = -[rl_bar ; rA];
    solution(p) = L'\(D\(L\rhs2(p)));

    dx = solution(1:length(x))';
    dy = solution(length(x)+1:length(x)+length(y))';
    
    dz = -(S\Z)*C'*dx+(S\Z)*(rC-Z\rSZ_bar);
    ds = -(Z\rSZ_bar)-(Z\(S*dz));
    
    % Update alpha
    dZS = [ dz ; ds ] ;
    alphas = (-[z; s]./dZS) ;
    alpha = min([1;alphas(dZS<0)]);

    % Update of position
    x = x+eta*alpha*dx;
    y = y+eta*alpha*dy;
    z = z+eta*alpha*dz;
    s = s+eta*alpha*ds;

    S = diag(s);
    Z = diag(z);
    e = ones(size(s));

    rL = H*x+g-A*y-C*z;
    rA = -A'*x+b;
    rC = -C'*x+s+d;
    rSZ = S*Z*e;

    dualGap = (z'*s)./m;

    % Converged
    Converged = (norm(rL) <= tol) && ...
            (norm(rA) <= tol) && ...
            (norm(rC) <= tol) && ...
            (abs(dualGap) <= tol) && ...
            (norm(s) <= tol);
%%

hist = vertcat(hist, 0.5*x'*H*x-g'*x);

end

hist
iter
x_sol = x;

end

