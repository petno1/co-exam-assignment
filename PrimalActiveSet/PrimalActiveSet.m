function [x_optimal] = PrimalActiveSet(H, g, A, b, C, dl, du, l, u)
%% Step 1 - prep the inputs
[n,~] = size(H);
% Convert dl <= C' x <= cu, l <= x <= u to c(x) = Cbar' x + dbar >= 0
Cbar = [C -C eye(n,n) -eye(n,n)];
dbar = [-dl; du; -l; u];
[~,m] = size(Cbar);
[~,j] = size(A);
% Tolerances
tol = 1.0e-8;

%% Step 2 - find x0 - some point which satisfies all constraints.
x0 = [A Cbar]' \ [b; dbar];

%% Step 3 - set up empty working/active sets

Wset = zeros(0,1);
IWset = [1:j+m]';
Wset = union(Wset, 1:j) ;  % add the integers 1-j to the set

lagrangeMultipliers = zeros(m+j,1);
mu =  zeros(j,1);
x = x0;

% QP data
gk = H*x + g;
nablaxL = gk - [A Cbar]*lagrangeMultipliers;
c = Cbar'*x + dbar;                   

% Check if the initial point is optimal
kktConditions = (norm(nablaxL,'inf') < tol); 

%% Main loop
maxit = 100;
it = 0;
    while ( ~kktConditions && (it < maxit) );
        it = it + 1
        % Solve equality constrained QP
        % join constraints
        AC = [A Cbar];
        bd = [b; dbar];

        Cw = AC(:,Wset);
        dw = zeros(size(Cw,2),1);
        [p,lagrangeMultipliers] = EqualityQPsubproblem(H,gk,Cw,dw);

        if(norm(p,'inf')>tol) % p is non-zero
            % find binding constraint (if any)
        alpha = 1.0;
        idc = -1;
        nIWset = size(IWset,1);
        for i = 1:nIWset
            pA = AC(:,IWset(i))'*p;
            if pA < 0.0
                alphapA = - c(IWset(i),1)/pA;
                if alphapA < alpha
                    alpha = alphapA;
                    idc = i;
                end
            end
        end
        % Take step, update data and working set
        x = x + alpha*p;
        gk = H*x + g;
        c = Cbar'*x + dbar;   
        if idc > 0
            Wset = [Wset; IWset(idc)];
            IWset = [IWset(1:idc-1); IWset(idc+1:end)];
        end
    else % p is zero
        % find minimum lambda
        idlambdaWset = -1;
        minlambdaWset = 0.0;
        nWset = size(Wset,1);
        for i=1:nWset
            if lagrangeMultipliers(i) < minlambdaWset 
                idlambdaWset = i;
                minlambdaWset = lagrangeMultipliers(i);
            end
        end
        if idlambdaWset > 0 % update the working set, x = x
            % if minimum lambda < 0 remove constraint from working set
            IWset = [IWset; Wset(idlambdaWset,1)];
            Wset = [Wset(1:idlambdaWset-1,1); Wset(idlambdaWset+1:end,1)];
        else % optimal solution found
            kktConditions = 1; %true
            xopt = x;
            lambdaopt = zeros(m,1);
            lagrangeMultipliers(Wset,1) = lagrangeMultipliers;
        end
    end
end

if ~kktConditions
    xopt = [];
    lambdaopt = [];
    Wset = [];
end

