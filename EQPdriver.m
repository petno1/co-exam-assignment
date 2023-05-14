clear;

[H, g, A, b] = GenerateRandomEQP(5,0.15);
[x_true, fval_true] = quadprog(H,g',[],[],A,b);

[x1] = EQPsolver(H,g,A,b,"LDLdense")
fval1 = calculateFval(x1, H, g)

[x2] = EQPsolver(H,g,A,b,"LDLsparse")
fval2 = calculateFval(x2, H, g)

[x3] = EQPsolver(H,g,A,b,"LUdense")
fval3 = calculateFval(x3, H, g)

[x4] = EQPsolver(H,g,A,b,"LUdense")
fval4 = calculateFval(x4, H, g)

[x5] = EQPsolver(H,g,A,b,"NullSpace")
fval5 = calculateFval(x5, H, g)

[x6] = EQPsolver(H,g,A,b,"RangeSpace")
fval6 = calculateFval(x6, H, g)

function fval = calculateFval(x, H, g)
    fval = 0.5 * x' * H * x + g' * x;
end