function [H, g, A, b, C, dl,du, l, u] = RandomQPGenerator(n, m, alpha, density)
    % Inputs:
    % n - integer, the size of the problem
    % m - integer, the number of constraints
    % alpha - scalar, a scaling factor for the identity matrix.
    % density - scalar, the density of the sparse matrices M, A and C. 
    % Outputs:
    % H - n x n positive definite matrix
    % g - n x 1 random vector
    % A - n x m sparse matrix
    % b - m x 1 random vector
    % C - n x m sparse matrix
    % du/dl - m x 1 random vector
    % l - n x 1 lower bounds
    % u - n x 1 upper bounds  
    
    M = sprand(n, n, density);
    H = M * M' + alpha * eye(n); 
    g = rand(n, 1); % Generate a random vector g
    x_rand = randn(n,1); % pick a random point
    
    % Generate random matrices A and C
    A = sprand(n, m, density);
    C = sprand(n, m, density);
    
    % Generate random vectors b dl du such that the random point will not be
    % out of bounds
    dl = randn(m,1);
    du = dl+C'*(x_rand);
    b = -A'*x_rand;
    
    % Generate random upper and lower bounds
    l = -ones(n,1)*2;
    u = ones(n,1)*2;
end
