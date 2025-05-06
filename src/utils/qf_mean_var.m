function [m,v] = qf_mean_var(mu,Cov,Q)

% [m,v] = f_mean_var(mu,Cov,LAM) returns mean m and variance v of the
% quadratic form y = x'*Q*x, where x are Gaussian random variables with
% mean vector mu and covariance matrix Cov
%
% Inputs:
%
%   mu: mean (column) vector of normal variables of quadratic form
%
%   Cov: covariance matrix of normal variables
%
%   Q: symmetric matrix of quadratic coefficients
%
% Ouputs:
%
%   m: expected value (mean) of quadratic form
%
%   v: variance of quadratic form
%
% Author: Paolo Manfredi
% Affiliation: Politecnico di Torino
% Date: February 2025

m = trace(Q*Cov) + mu'*Q*mu;

if nargin>1
    v = 2*trace(Q*Cov*Q*Cov) + 4*mu'*Q*Cov*Q*mu;
end