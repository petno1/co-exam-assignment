function [x_sol] = InteriorPointQP(H, g, A, b, C, d, x0)

% Convert Cx >= d into an equality constraint using a slack variable
%   A' x = b
%   C' x + d + s = 0
%   s > 0

% Sets constants for the algorithm
m = length(d);
tol = 0.00000001;
eta = 0.995;
x0 = x0 +tol;
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
rC = -C'*x+s-d;
rSZ = S*Z*e;
h_bar = H + C*(S\Z)*C';
rl_bar =  rL-C*(S\Z)*(rC-Z\rSZ);
KKT = [h_bar -A; -A' zeros(size(A,2))];
[L,D,p] = ldl(KKT,"lower",'vector');    
rhs = -[rl_bar;rA];
solution(p) = L'\(D\(L\rhs(p)));

%% Affine Direction
dx_affine = solution(1:length(x))';
dz_affine =-(S\Z)*C'*dx_affine+(S\Z)*(rC-Z\rSZ);
ds_affine = (-Z)\rSZ-(Z\S)*dz_affine;
z = max(1, abs(z0 + dz_affine));
s = max(1, abs(s0 + ds_affine));

dualGap0 = (z'*s)/(m);
dualGap = (z'*s)/(m);

% Converged
Converged = (norm(rL,inf) <= tol) && ...
            (norm(rA,inf) <= tol) && ...
            (norm(rC,inf) <= tol) && ...
            (abs(s) <= tol) && ...
            (max(abs(s.*z)) <= tol);
%%       
max_iter = 20;
iter = 0;
hist = [];

while ~Converged && (iter<max_iter)

    iter = iter+1;

    h_bar = H + C*(S\Z)*C';
    rl_bar =  rL-C*(S\Z)*(rC-Z\rSZ);
    KKT = [h_bar -A; -A' zeros(size(A,2))];
    [L,D,p] = ldl(KKT,"lower",'vector');    
    rhs = -[rl_bar;rA];
    solution(p) = L'\(D\(L\rhs(p)));

    %% Affine Direction
    dx_affine = solution(1:length(x))';
    dz_affine =-(S\Z)*C'*dx_affine+(S\Z)*(rC-Z\rSZ);
    ds_affine = (-Z)\rSZ-(Z\S)*dz_affine;

   % Calculate alpha
    dZS = [ dz_affine ; ds_affine ] ;
    alphas = (-[z; s]./dZS) ;
    affineAlpha = min([1;alphas(dZS<0)]);

    %% Duality gap and centering parameter
    affineDualGap = (((z+affineAlpha*dz_affine)'*(s+affineAlpha*ds_affine)))*(1/m);
    sigma = (affineDualGap./dualGap)^3;

    %% Affine-Centering-Correction Direction
    rSZ_bar = rSZ + (ds_affine.*dz_affine.*e) - (sigma.*affineDualGap*e);
    rl_bar = rL- C*(S\Z)*(rC-Z\rSZ_bar);


    rhs2 = -[rl_bar ; rA];
    solution(p) = L'\(D\(L\rhs2(p)));

    dx = solution(1:length(x))';
    dy = solution(length(x)+1:length(solution))';
    
    dz = -(S\Z)*C'*dx+(S\Z)*(rC-Z\rSZ_bar);
    ds = (-Z)\rSZ_bar-Z\S*dz;
    
    % Update alpha
    dZS = [dz;ds] ;
    alphas = (-[z;s]./dZS) ;
    alpha = min([1;alphas(dZS<0)]);

    alpha_bar = alpha*eta;
    % Update of position
    x = x+alpha_bar*dx;
    y = y+alpha_bar*dy;
    z = z+alpha_bar*dz;
    s = s+alpha_bar*ds;

    S = diag(s);
    Z = diag(z);
    e = ones(size(s));

    rL = H*x+g-A*y-C*z;
    rA = b - A'*x;
    rC = -C'*x+s-d;
    rSZ = S*Z*e;
    
    dualGap = (z'*s)/(m);

    Converged = (norm(rL,inf) <= tol*max(1, max(union(union(union(A,g),H),C)))) && ...
            (norm(rA,inf) <= tol*max(1, max(union(A', b)))) && ...
            (norm(rC,inf) <= tol*max(1, max(union(d, C')))) && ...
            (dualGap <= tol*0.01*dualGap0);


hist = vertcat(hist, 0.5*x'*H*x-g'*x);

end

hist
x_sol = x;
disp(iter)
if iter == max_iter
    disp("hej")
end
end
