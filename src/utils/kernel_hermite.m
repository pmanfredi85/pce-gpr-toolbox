function Omega = kernel_hermite(x,xp,rho)

% Univariate Hermite (Mehler) kernel
% 
% Omega = kernel_hermite(x,xp,rho) accepts an Mx1 input x and an Nx1 input
% xp, and returns the MxN kernel matrix. The hyperparameter rho is a scalar
% in the interval (0,1).
%
% Author: Paolo Manfredi
% Affiliation: Politecnico di Torino
% Date: February 2025

Omega = 1/sqrt(1-rho^2)*exp(-(rho^2*(x.^2+xp'.^2)-2*rho*x.*xp')/(2*(1-rho^2)));