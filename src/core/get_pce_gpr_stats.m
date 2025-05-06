function [mu,sigma2,St,var_mu,var_sigma2,var_St,pdf_mu,pdf_sigma2,pdf_St] = get_pce_gpr_stats(M_pce_gpr)

% [mu,sigma2,St,var_mu,var_sigma2,var_St,pdf_mu,pdf_sigma2,pdf_St] =
% get_pce_gpr_stats(M_pce_gpr) computes PCE-based statistical information
% from the PCE-GPR model M_pce_gpr
%
%   Note: this function requires the Generalized chi-square distribution
%   toolbox by Abhranil Das - version 2.3.0 to handle the distribution of
%   output variance and Sobol' indices
%
%   https://it.mathworks.com/matlabcentral/fileexchange/85028-generalized-chi-square-distribution
%
%   A. Das, "New methods for computing the generalized chi-square
%   distribution." arXiv preprint arXiv:2404.05062 (2024).
% 
% Inputs:
%
%   M_pce_gpr: a PCE-GPR model structure trained with train_pce_gpr and
%   including PCE coefficients obtained with get_pce_gpr_coeff
%
% Outputs:
%
%   mu: output mean
%
%   sigma2: output variance
%
%   St: Sobol' indices (non-normalized); a vector of length d, where d is
%   the number of input parameters
%
%   var_mu: variance of the output mean
%
%   var_sigma2: variance of the output variance
%
%   var_St: variances of the Sobol' indices; a vector of length d
%
%   pdf_mu: distribution of the output mean; a NormalDistribution object
%
%   pdf_sigma2: distribution of the output variance;
%   GeneralizedChiSquareDistribution object
%
%   pdf_St: distributions of the Sobol' indices; an array of length d of
%   GeneralizedChiSquareDistribution objects
%
% Author: Paolo Manfredi
% Affiliation: Politecnico di Torino
% Date: March 2025

% number of input parameters
d = size(M_pce_gpr.x_train,2);


%% output mean

% expected value
mu = M_pce_gpr.pce_coeff(1);

% variance of output mean
if nargout>3
    var_mu = M_pce_gpr.cov_pce_coeff(1,1);
end

% pdf of output mean
if nargout>6
    pdf_mu = makedist("Normal","mu",mu,"sigma",sqrt(var_mu));
end

%% output variance

if nargout>1
    % quadratic form
    mu_qf = M_pce_gpr.pce_coeff(2:end);
    Cov_qf = M_pce_gpr.cov_pce_coeff(2:end,2:end);
    Q_qf = eye(size(Cov_qf));

    % expected value & variance
    [sigma2,var_sigma2] = qf_mean_var(mu_qf,Cov_qf,Q_qf);
end

% pdf
if nargout>7
    [w,k,lambda,m,s] = norm_quad_to_gx2_params(mu_qf,Cov_qf,struct('q2',Q_qf,'q1',zeros(length(mu_qf),1),'q0',0));
    pdf_sigma2 = GeneralizedChiSquareDistribution(real(w),real(k),real(lambda),real(m),real(s));
end

if nargout>2

    St = zeros(d,1);
    var_St = zeros(d,1);

    % Sobol' indices
    for jj = 1:d
        index = M_pce_gpr.kvec(:,jj)>0;

        % quadratic form
        mu_qf = M_pce_gpr.pce_coeff(index);
        Cov_qf = M_pce_gpr.cov_pce_coeff(index,index);
        Q_qf = eye(size(Cov_qf));

        % expected value & variance
        [St(jj),var_St(jj)] = qf_mean_var(mu_qf,Cov_qf,Q_qf);

        % pdf
        if nargout>8
            [w,k,lambda,m,s] = norm_quad_to_gx2_params(mu_qf,Cov_qf,struct('q2',Q_qf,'q1',zeros(length(mu_qf),1),'q0',0));
            pdf_St(jj) = GeneralizedChiSquareDistribution(real(w),real(k),real(lambda),real(m),real(s));
        end
    end
end