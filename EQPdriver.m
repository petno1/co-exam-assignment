[H, g, A, b] = GenerateRandomEQP(5,0.15);

[x1] = EQPsolver(H,g,A,b,"LDLdense");
[x2] = EQPsolver(H,g,A,b,"LDLsparse");

fval1 = 0.5*x1'*H*x1+g'*x1
fval2 = 0.5*x2'*H*x2+g'*x2

[~, fval_true] = quadprog(H,g',[],[],A,b)
