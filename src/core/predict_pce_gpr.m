function [y_gpr,y_pce,cov_gpr,cov_pce] = predict_pce_gpr(M_pce_gpr,x)

% [y_gpr,y_pce,Cov_gpr,Cov_pce] = predict_pce_gpr(M_pce_gpr,x) evaluates
% the PCE-GPR model M_pce_gpr at points x
%
% Inputs:
%
%   M_pce_gpr: structure describing the PCE-GPR model, returned by
%   train_pce_gpr
%
%   x: matrix of points at which to evaluate the model of size n x d, where
%   d is the number of dimensions
%
% Outputs:
%
%   y_gpr: prediction of the GPR model (dual space kernel formulation)
%
%   y_pce: prediction of the PCE model (truncated primal space formulation)
%
%   cov_gpr: posterior covariance of predictions from the GPR model
%
%   cov_pce: posterior covariance of predictions from the PCE model
%
%   Note: The PCE describes only a projection onto a finite set of basis
%   functions, therefore y_pce has lower accuracy than y_gpr. It is
%   calculated only if the PCE coefficients are available in the model
%   structure (they need to be computed on beforehand using
%   get_pce_gpr_coeff). The same consideration applies to cov_gpr and
%   cov_pce
%
% Author: Paolo Manfredi
% Affiliation: Politecnico di Torino
% Date: March 2025


%% reading necessary input parameters

x_train = M_pce_gpr.x_train;
kernel = M_pce_gpr.kernel;
alp = M_pce_gpr.alpha;
distr = M_pce_gpr.distr;
basis_functions = M_pce_gpr.basisfun;

[~,d1] = size(x_train);
[N,d2] = size(x);

if d2~=d1
    error('Number of dimensions of x and x_train do not match');
end

% standardize inputs
x = standardize_input(x,distr);


%% evaluate model

% GPR model
Omega_val = kernel(x,x_train);

y_gpr = Omega_val*alp;

% PCE model
if nargout>1
    if isfield(M_pce_gpr,'pce_coeff')
        kvec = M_pce_gpr.kvec;
        c = M_pce_gpr.pce_coeff;

        Psi = evalBasisFunctions_heterogeneous(basis_functions,kvec,x);
        y_pce = Psi*c;
    else
        y_pce = [];
    end
end

% posterior covariance
if nargout>2
    tau = M_pce_gpr.tau;
    L = M_pce_gpr.L;
    sigma2_tot = M_pce_gpr.sigma^2 + M_pce_gpr.sigma_n^2;

    rtilde = (1-tau)*kernel(x_train,x);
    invLr = L\rtilde;

    cov_gpr = sigma2_tot*(kernel(x,x) - invLr'*invLr);

    if nargout>3 && isfield(M_pce_gpr,'pce_coeff')
        cov_pce = zeros(N);

        % Bienaymé's identity
        for ii = 1:N
            for jj = 1:N
                cov_pce(ii,jj) = sum(sum(Psi(ii,:).*M_pce_gpr.cov_pce_coeff.*Psi(jj,:)'));
            end
        end
    else
        cov_pce = [];
    end
end
