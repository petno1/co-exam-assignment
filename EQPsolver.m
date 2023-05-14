function [x, lambda] = EQPsolver(H, g, A, b, type)
    n = size(H, 1);
   
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
        [K, d] = sparseConstructKKTSystem(H, g, A, b);
        [L, U, p] = lu(K, 'vector'); 
        solution(p) = U\(L\( d(p)));
        x = solution(1:n)';
        lambda = solution(n+1:end)';
    end 
    
end
