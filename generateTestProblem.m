function [H, g, A, b] = generateTestProblem(n)

    % generates a random quadratic programming problem with n variables and m equality constraints
    % input n: The number of variables in the problem

    m = randi([1, n-1]); % Generate a random integer between 1 and n-1

    % Generate a random positive definite matrix H
    H = randn(n);
    H = 0.5 * (H + H'); % Make H symmetric
    H = H + eye(n)*n; % Ensure H is positive definite

    % Generate a random matrix A
    A = randn(n, m);
    A = -A; % Ensure A*x <= 0

    b = randi(n)*ones(m, 1);
    g = randi(n)*ones(n, 1);
  
end