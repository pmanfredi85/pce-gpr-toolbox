function Omega = kernel_legendre(x,xp,rho)

% Univariate Legendre kernel
% 
% Omega = kernel_legendre(x,xp,rho) accepts an Mx1 input x and an Nx1 input
% xp, and returns the MxN kernel matrix. The hyperparameter rho is a scalar
% in the interval (0,1).
%
% Author: Paolo Manfredi
% Affiliation: Politecnico di Torino
% Date: February 2025

K_fun = @(k) ellipke(k.^2);
a_fun = @(x,xp,rho) 1 - 2*x.*xp'*rho + rho^2;
b_fun = @(x,xp,rho) -2*rho*sqrt(1-x.^2).*sqrt(1-xp'.^2);
Omega = 2*K_fun(sqrt(2*b_fun(x,xp,rho)./(b_fun(x,xp,rho)-a_fun(x,xp,rho))))./(pi*sqrt(a_fun(x,xp,rho)-b_fun(x,xp,rho)));