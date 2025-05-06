function Omega = kernel_heterogeneous(x,xp,rho,distr)

% Omega = kernel_heterogeneous(x,xp,rho,distr) evaluates, at points
% specified by x and xp, the separable kernel constructed as the product of
% univariate kernels associated to the probability distributions specified
% by distr
%
% Inputs:
%
%   x: matrix of first input samples of size M x d, where M is the number
%   of samples and d is the number of dimensions
%
%   xp: matrix of second input samples of size N x d, where N is the number
%   of samples and d is the number of dimensions
%
%   rho: kernel hyperparameters. It must be a scalar for an isotropic
%   kernel, or a vector of length d for an anisotropic kernel
%
%   distr: structure array of length d specifying the distributions and
%   hence the kernel functions to be used for each input dimension
%
% Output:
%
%   Omega: kernel matrix of size M x N, with the kernel evaluated at all
%   combinations of inputs x and xp
%
% Author: Paolo Manfredi
% Affiliation: Politecnico di Torino
% Date: February 2025

% processing inputs
[M,d1] = size(x);
[N,d2] = size(xp);

if d1~=d2
    error('x and xp appear to have different number of dimensions');
else
    d = d1;
end

% isotropic vs anisotropic kernel
if isscalar(rho)
    rho = rho*ones(1,d);
else
    if length(rho)~=d
        error('rho has incorrect size');
    end
end

% evaluate kernel
Omega = zeros(M,N,d);
for jj = 1:d
    switch lower(distr(jj).Type)
        case 'norm'
            Omega(:,:,jj) = kernel_hermite(x(:,jj),xp(:,jj),rho(jj));
        case 'unif'
            Omega(:,:,jj) = kernel_legendre(x(:,jj),xp(:,jj),rho(jj));
        otherwise
            error('Wrong distribution type');
    end
end

Omega = prod(Omega,3);