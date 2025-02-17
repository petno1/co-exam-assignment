function [x, lambda] = EQPsolver(H, g, A, b, type)
    [~,n] = size(H);
    if strcmp(type, 'LDLdense')
        [K, d] = constructKKTSystem(H,g,A,b);
        [L,D,p] = ldl(K,'vector');
        solution(p) = L'\(D\(L\d(p)));
        x = solution(1:n)';
        lambda = solution(n+1:end)';
    end

     if strcmp(type, 'LDLsparse')
        [K, d] = sparseConstructKKTSystem(H, g, A, b);
        [L,D,p] = ldl(K,'vector');
        solution(p) = L'\(D\(L\d(p)));
        x = solution(1:n)';
        lambda = solution(n+1:end)';
    end
     
    if strcmp(type, 'LUdense')
        [K, d] = constructKKTSystem(H, g, A, b);
        [L, U, p] = lu(K, 'vector'); 
        solution(p) = U\(L\( d(p)));
        x = solution(1:n)';
        lambda = solution(n+1:end)';
    end 

     if strcmp(type, 'LUsparse')
        [K, d] = sparseConstructKKTSystem(H, g, A, b);
        [L, U, p] = lu(K, 'vector'); 
        solution(p) = U\(L\( d(p)));
        x = solution(1:n)';
        lambda = solution(n+1:end)';
     end 

    if strcmp(type, "RangeSpace")
        R=chol(H);
        mu=R'\g ;
        Hg=R\mu;
        mu=R'\A';
        HA=R\mu;   
        lambda = (A*HA) \(b+A*Hg) ;
        x = HA*lambda-Hg;
    end

    if strcmp(type, "NullSpace")
        [Q,R_bar] = qr(A');
        m = size(R_bar,2);
        Q1 = Q(:,1:m);
        Q2 = Q(:,m+1:end);
        R = R_bar(1:m,1:m);
        xy = R'\b;
        xz = (Q2'*H*Q2)\(-Q2'*(H*Q1*xy+g));
        x = Q1*xy+Q2*xz;
        lambda = R\(Q1'*(H*x+g));
    end

end
