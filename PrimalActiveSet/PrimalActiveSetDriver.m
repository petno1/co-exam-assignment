[x] = PrimalActiveSet(H, g, A, b, C, dl, du, l, u)
[x2] = quadprog(H, g', C',du-dl, A', b, l, u)