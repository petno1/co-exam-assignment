function [x_sol] = InteriorPointQP(H, g, A, b, C, d, x0)

% Convert Cx >= d into an equality constraint using a slack variable
%   A' x = b
%   C' x + d + s = 0
%   s > 0

% Sets constants for the algorithm
m = length(d);
tol = 0.000001;
eta = 0.995;

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

dualGap = (z'*s)/(m);

% Converged
Converged = (norm(rL,inf) <= tol) && ...
            (norm(rA,inf) <= tol) && ...
            (norm(rC,inf) <= tol) && ...
            (abs(s) <= tol);
%%       
max_iter = 200;
iter = 0;
hist = [];

while ~Converged && (iter<max_iter)

    iter = iter+1;

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

   % Calculate alpha
    affineDualGap = 1;
    idelta_x_z = find(dz_affine<0);

    if (isempty(idelta_x_z) == 0)
        affineDualGap = min(affineDualGap, min(-z(idelta_x_z)./dz_affine(idelta_x_z)));
    end
    idelta_x_s = find(ds_affine<0);
    if (isempty(idelta_x_s) == 0)
        affineDualGap = min(affineDualGap, min(-s(idelta_x_s)./ds_affine(idelta_x_s)));
    end

    %% Duality gap and centering parameter
    affineDualGap = ((z+affineDualGap*dz_affine)'*(s+affineDualGap*ds_affine))/(m);
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
    
    % Update alpha
    alpha = 1;
    idelta_x_z = find(dz<0);
    
    if (isempty(idelta_x_z) == 0)
        alpha = min(alpha, min(-z(idelta_x_z)./dz(idelta_x_z)));
    end

    idelta_x_s = find(ds<0);

    if (isempty(idelta_x_s) == 0)
        alpha = min(alpha, min(-s(idelta_x_s)./ds(idelta_x_s)));
    end

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

    dualGap = (z'*s)/(m);

    Converged = (norm(rL,inf) <= tol) && ...
            (norm(rA,inf) <= tol) && ...
            (norm(rC,inf) <= tol) && ...
            (abs(s) <= tol);

hist = vertcat(hist, 0.5*x'*H*x-g'*x);

end

hist
x_sol = x;

end

