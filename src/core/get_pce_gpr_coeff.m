function [M_pce_gpr,c] = get_pce_gpr_coeff(M_pce_gpr,kvec)

% [M_pce_gpr,c] = get_pce_gpr_coeff(M_pce_gpr,kvec) retrieves the
% polynomial chaos expansion (PCE) coefficients of the basis functions
% indicated by kvec from the PCE-GPR model M_pce_gpr
%
% Inputs:
%
%   M_pce_gpr: PCE-GPR model structure trained using train_pce_gpr
%
%   kvec: array of multi-indices defining the PCE basis functions
%
% Ouputs:
%
%   M_pce_gpr: updated PCE-GPR model structure including the following
%   fields:
%
%       M_pce_gpr.kvec: array of PCE multi-indices (input)
%
%       M_pce_gpr.pce_coeff: array of PCE coefficients corresponding to kvec
%
%       M_pce_gpr.corr_pce_coeff: correlation matrix of the PCE coefficients
%
%       M_pce_gpr.cov_pce_coeff: covariance matrix of the PCE coefficients
%
%   c: vector of PCE coefficients (same as M_pce_gpr.pce_coeff)
%
% Author: Paolo Manfredi
% Affiliation: Politecnico di Torino
% Date: March 2025


% retrieve relevant information for internal use
[~,d] = size(M_pce_gpr.x_train); % number of input parameters

% vector of rescale factors
lambda_vec = zeros(size(kvec));

if isscalar(M_pce_gpr.rho)
    rho_vec = M_pce_gpr.rho*ones(1,d);
else
    rho_vec = M_pce_gpr.rho;
end

for jj = 1:d
    switch lower(M_pce_gpr.distr(jj).Type)
        case 'norm'
            lambda_vec(:,jj) = rho_vec(jj).^kvec(:,jj);

        case 'unif'
            lambda_vec(:,jj) = rho_vec(jj).^kvec(:,jj)./(2*kvec(:,jj)+1);
    end
end

lambda_vec = prod(lambda_vec,2);

% PCE basis functions at training samples
Psi = evalBasisFunctions_heterogeneous(M_pce_gpr.basisfun,kvec,M_pce_gpr.x_train);

% PCE coefficients (dual space to primal space conversion)
c = (Psi'.*lambda_vec)*M_pce_gpr.alpha;

% correlation matrix of PCE coefficients
LinvPsiLam = M_pce_gpr.L\(Psi.*lambda_vec');
M_pce_gpr.corr_pce_coeff = (1-M_pce_gpr.tau)*(diag(lambda_vec) - (1-M_pce_gpr.tau)*(LinvPsiLam'*LinvPsiLam));

M_pce_gpr.kvec = kvec;
M_pce_gpr.pce_coeff = c;

% covariance matrix of PCE coefficients
M_pce_gpr.cov_pce_coeff = (M_pce_gpr.sigma^2+M_pce_gpr.sigma_n^2)*M_pce_gpr.corr_pce_coeff;