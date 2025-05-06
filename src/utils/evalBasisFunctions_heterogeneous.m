function Psi = evalBasisFunctions_heterogeneous(basisfun,kvec,x)

% Psi = evalBasisFunctions_heterogeneous(basisfun,kvec,x) evaluates
% the multivariate basis functions specified by basisfun (one-dimensional
% basis) and multi-index matrix kvec at points x
%
% Inputs:
%
%   basisfun(k,x): function handle to one-dimensional basis functions.
%   basisfun must accept a vector of points x and an integer order k
%
%   kvec: matrix of multi-indices of multivariate basis functions. The
%   multivariate basis functions are constructed as the product of
%   one-dimensional ones. Each row indicates the degree of the
%   corresponding function in each dimension
%
%   x: matrix of points at which to evaluate the basis functions. Note: the
%   number of columns must equal the number of dimensions and hence the
%   number of columns of kvec
%
% Output:
%
%   Psi: matrix of the basis functions evaluated at points x. Psi is a L x
%   K matrix, where L is the number of points (i.e., rows of x) and K is
%   the number of basis functions (i.e., rows of kvec)
%
% Author: Paolo Manfredi
% Affiliation: Politecnico di Torino
% Date: February 2025

[K,d] = size(kvec);
[L,d2] = size(x);

if d~=d2
    error('kvec and x have incompatible number of columns');
end

% Building matrix Psi
Psi = ones(L,K);
for k = 1:K
    for jj = 1:d
        Psi(:,k) = Psi(:,k) .* basisfun{jj}(kvec(k,jj),x(:,jj));
    end
end