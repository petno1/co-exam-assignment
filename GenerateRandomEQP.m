function [H, g, A, b] = GenerateRandomEQP(n, density)

   m = round(rand*n);
    
   M = sprandn(n,n,density);
   H = M*M' + 0.001*eye(n); 
    
   A = sprandn(m,n,density);
   A(1:m+1:end) = 1;

   x_random = randn(n,1);
   lambda_random = randn(m,1);

   KKT = [H, -A'; -A, zeros(m, m)];
   solution = -KKT * [x_random; lambda_random];
   g = -solution (1 : n);
   size(g)
   b = -solution (n+1:end);

end 


