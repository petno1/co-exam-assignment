function [xopt,lagrangeopt,Wset,it] = PrimalActiveSet(H, g, A, b, Cbar, dbar, x0)

%% Step 1 - prep the inputs
[n,~] = size(H);
[~,m] = size(Cbar);
[~,j] = size(A);
tol = 1.0e-10;
x = x0;

%% Step 2 - set up empty working/active sets
Wset = zeros(0,1);
IWset = [1:j+m]';
Wset = union(Wset, 1:j);  % add the integers 1-j to the set
lagrangeMultipliers = zeros(m+j,1);



% QP data
gk = H*x + g;
nablaxL = gk - [A Cbar]*lagrangeMultipliers;
c = [A Cbar]'*x + [b; dbar];                   

% Check if the initial point is optimal
kktConditions = (norm(nablaxL,'inf') < tol); 

%% Main loop
maxit = 100;
it = 0;
while (~kktConditions && (it < maxit))      
    Wset
    it = it + 1;
    % Solve equality constrained QP
    % join constraints
    AC = [A Cbar];
    Cw = AC(:,Wset);
    dw = zeros(size(Cw,2),1);
    [p,lagrangeMultipliers] = EqualityQPsubproblem(H,gk,Cw,dw);
    if(norm(p,'inf')>tol) % p is non-zero
        % find binding constraint (if any)
        alpha = 1;
        idc = -1;
        nIWset = size(IWset,1);
        for i = j+1:nIWset
            pA = AC(:,IWset(i))'*p;
            if pA < 0.0
                alphapA = -c(IWset(i),1)/pA;
                if alphapA < alpha
                    alpha = alphapA;
                    idc = i;
                end
            end
        end
        % Take step, update data and working set
        x = x + alpha*p;
        gk = H*x + g;
        c = [A Cbar]'*x + [b; dbar];                   
        if idc > 0
            Wset = [Wset; IWset(idc)];
            IWset = [IWset(1:idc-1); IWset(idc+1:end)];
        end
    else % p is zero
        % find minimum lambda
        idLagrangeWset = -1;
        minLagrangeWset = 0.0;
        nWset = size(Wset,1);
        for i=j+1:nWset
            if lagrangeMultipliers(i) < minLagrangeWset 
                idLagrangeWset = i;
                minLagrangeWset = lagrangeMultipliers(i);
            end
        end
        if idLagrangeWset > 0 % update the working set, x = x
            % if minimum lambda < 0 remove constraint from working set
            IWset = [IWset; Wset(idLagrangeWset,1)];
            Wset = [Wset(1:idLagrangeWset-1,1); Wset(idLagrangeWset+1:end,1)];
        else % optimal solution found
            disp('BIG BONG DONE')
            kktConditions = 1; %true
            xopt = x;
            lagrangeopt = zeros(m,1);
            lagrangeopt(Wset,1) = lagrangeMultipliers;
        end
    end
end

if ~kktConditions
    xopt = [];
    lagrangeopt = [];
    Wset = [];
end

