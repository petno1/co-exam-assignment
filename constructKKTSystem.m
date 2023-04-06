function [K,d] = constructKKTSystem(H,g,A,b,varargin)
%   Inputs:
%   H: n-by-n symmetric matrix of quadratic coefficients
%   g: n-by-1 column vector of linear coefficients
%   A: m-by-n matrix of linear equality constraints
%   b: m-by-1 column vector of equality constraint values
%   sparse: optional flag to indicate whether to use sparse matrices (default: false)

% Parse optional parameter
use_sparse = false;
if ~isempty(varargin)
    use_sparse = varargin{1};
end

% Determine matrix size
[m,m] = size(A);

% Construct KKT system
if use_sparse
    K = sparse([H, -A; -A', zeros(m,m)]);
else
    K = [H, -A; -A', zeros(m,m)];
end
d = [-g; -b];
end