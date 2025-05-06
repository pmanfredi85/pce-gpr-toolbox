function M_pce_gpr = train_pce_gpr(x_train,y_train,distr,varargin)

% M_pce_gpr = train_pce_gpr(x_train,y_train,distr,varargin) uses UQLab to
% train a PCE-GPR model of an N-D function based on training data
% (x_train,y_train) using a Gaussian process (GPR) forumulation.
%
% Inputs:
%
%   x_train: matrix of input training samples of size n x d, where n is the
%   number of samples and d is the number of dimensions
%
%   y_train: vector of output observations at the training samples, of size
%   n x 1, where n is the number of samples
%
%   distr: structure array describing the probability distributions of
%   input parameters. The structure has two fields: Type and Parameters. A
%   different distribution for each input parameter can be specified. At
%   the moment, only Gaussian and uniform distributions are supported
%
%       Type: distribution type
%           'norm' = Gaussian
%           'unif' = uniform
%
%       Parameters: parameters of distribution
%           If Type = 'norm', Parameters must be an array with mean and
%           standard deviation;
%           If Type = 'unif', Parameters must be an array with lower and
%           upper bound
%
% Optional Name-Value arguments:
%
%   IsotropicKernel: specifies whether to use an isotropic or anisotropic
%   kernel
%
%       false (default): use an anisotropic kernel (different rho for each input)
%
%       true: use an isotropic kernel (same rho for all inputs)
%
%   Noise: specifies whether to include noise on training data
%
%       false (default): neglect noise (sets noise variance to zero)
%
%       true: do consider noise
%
%   EstimationMethod: hyperparameter estimation method
%
%       'ML' (default): maximum likelihood estimation
%
%       'CV': leave-one-out cross-validation
%
%   UQLabOpts: structure with advanced UQLab options for hyperparameter
%   estimation. Allows specifying options below (for additional details,
%   please refer to UQLab manual):
%
%       .CV.LeaveKOut: left-out samples in cross-validation. Default: 1 (UQLab default, = leave-one-out)
%       .Optim.InitialValue: intial value(s) of rho. Default: 0.1 (scalar or vector, depending on IsotropicKernel)
%       .Optim.Bounds: bounds for rho. Default: [1e-6,1] (scalar or vector, depending on IsotropicKernel). Note that rho must be bounded between 0 and 1.
%       .Optim.Display: verbosity of optimization process. Default: 'iter'
%       .Optim.MaxIter: maximum number of iterations. Default: 20 (UQLab default)
%       .Optim.Tol: convergence tolerance. Default: 1e-4 (UQLab default)
%       .Optim.Method: optimization method. Default: 'HCMAES'
%       .Optim.BFGS: BFGS specific options
%       .Optim.GA: GA specific options
%       .Optim.HGA: HGA specific options
%       .Optim.CMAES: CMAES specific options
%       .Optim.HCMAES: HCMAES specific options
%
% Outputs:
%
%   M_pce_gpr: structure with information concerning the trained model
%   M_pce_gpr.x_train: standardized training samples (rescaled input)
%   M_pce_gpr.y_train: training observations (input)
%   M_pce_gpr.distr: structure array with input distributions (input)
%   M_pce_gpr.basisfun: cell array with basis functions for each input dimension
%   M_pce_gpr.kernel: kernel function
%   M_pce_gpr.L: Cholesky factor of matrix Rtilde (correlation matrix of the training samples)
%   M_pce_gpr.rho: kernel hyperparameters
%   M_pce_gpr.tau: noise ratio
%   M_pce_gpr.sigma: kernel standard deviation
%   M_pce_gpr.sigma_n:noise standard deviation
%   M_pce_gpr.alpha: GPR coefficients
%   M_pce_gpr.model: UQLab model object
%
% Author: Paolo Manfredi
% Affiliation: Politecnico di Torino
% Date: February 2025


%% defining parameters

% optimization settings
LogRhoLB = -6; % lower bound for logarithm of kernel hyperparameter(s)
LogRhoUB = 0; % upper bound for logarithm of kernel hyperparameter(s)
nugget = 1e-12; % nugget for covariance matrix


%% parse input

expectedDistrType = {'norm','unif'};
validDistr = @(x) isstruct(x);

defaultIsotropicKernel = false;
validIsotropicKernel = @(x) islogical(x);

defaultNoise = false;
validNoise = @(x) islogical(x);

defaultEstimationMethod = 'ML';
expectedEstimationMethod = {'CV','ML'};
validEstimationMethod = @(x) any(validatestring(x,expectedEstimationMethod));

defaultUQLabOpts = [];

ip = inputParser;
addRequired(ip,'x_train');
addRequired(ip,'y_train');
addRequired(ip,'distr',validDistr);
addParameter(ip,'IsotropicKernel',defaultIsotropicKernel,validIsotropicKernel);
addParameter(ip,'Noise',defaultNoise,validNoise);
addParameter(ip,'EstimationMethod',defaultEstimationMethod,validEstimationMethod);
addParameter(ip,'UQLabOpts',defaultUQLabOpts);

parse(ip,x_train,y_train,distr,varargin{:});

% read parameters
IsotropicKernel = ip.Results.IsotropicKernel;
Noise = ip.Results.Noise;
EstimationMethod = ip.Results.EstimationMethod;
UQLabOpts = ip.Results.UQLabOpts;

% further checks on inputs
[Ntrain,d] = size(x_train);

if length(y_train)~=Ntrain
    error('x_train and y_train are of incompatible size');
end

% check on distribution
if length(distr)~=d
    error('distr has incorrect size');
else
    for jj = 1:d
        any(validatestring(distr(jj).Type,expectedDistrType));
    end
end

% rescale input values to standardized distributions
x_train = standardize_input(x_train,distr);


%% Define kernel and basis functions

% define multivariate kernel
kernel = @(x,y,rho) kernel_heterogeneous(x,y,rho,distr);

% define basis functions
basis_functions = cell(d,1);
for jj = 1:d
    switch lower(distr(jj).Type)
        case 'norm'
            basis_functions{jj} = @(k,x) polyval(orthonormal_hermite(k),x);

        case 'unif'
            basis_functions{jj} = @(k,x) polyval(orthonormal_legendre(k),x);

    end
end


%% GPR model training

% uqlab initialization
uqlab
MetaOpts.Type = 'Metamodel';
MetaOpts.MetaType = 'Kriging';

% experimental design
MetaOpts.ExpDesign.X = x_train;
MetaOpts.ExpDesign.Y = y_train;

% trend (forced to be zero)
MetaOpts.Trend.Type = 'simple';
MetaOpts.Trend.CustomF = @(x) 0;

% kernel function (custom)
MetaOpts.Corr.Nugget = nugget;
MetaOpts.Corr.Handle = @(x1,x2,theta,options) kernel(x1,x2,10.^theta);

% estimation method (defined by input)
MetaOpts.EstimMethod = EstimationMethod;

% CV options
if isfield(UQLabOpts,'CV')
    MetaOpts.CV = UQLabOpts.CV;
end

% pass possible UQLab-specific options
if ~isfield(UQLabOpts,'Optim')
    UQLabOpts.Optim = []; % create a dummy .Optim field to ease subsequent checks
end
if isfield(UQLabOpts.Optim,'InitialValue')
    MetaOpts.Optim.InitialValue = log10(UQLabOpts.Optim.InitialValue);
else
    if IsotropicKernel
        MetaOpts.Optim.InitialValue = log10(0.1);
    else
        MetaOpts.Optim.InitialValue = log10(0.1)*ones(d,1);
    end
end
if isfield(UQLabOpts.Optim,'Bounds')
    if ~IsotropicKernel && size(UQLabOpts.Optim.Bounds,2)==1
        MetaOpts.Optim.Bounds = repmat(log10(UQLabOpts.Optim.Bounds),[1,d]);
    else
        MetaOpts.Optim.Bounds = log10(UQLabOpts.Optim.Bounds);
    end
else
    if IsotropicKernel
        MetaOpts.Optim.Bounds = [LogRhoLB,LogRhoUB]';
    else
        MetaOpts.Optim.Bounds = repmat([LogRhoLB,LogRhoUB]',[1,d]);
    end
end
if isfield(UQLabOpts.Optim,'Display')
    MetaOpts.Optim.Display = UQLabOpts.Optim.Display;
else
    MetaOpts.Optim.Display = 'iter';
end
if isfield(UQLabOpts.Optim,'MaxIter')
    MetaOpts.Optim.MaxIter = UQLabOpts.Optim.MaxIter;
end
if isfield(UQLabOpts.Optim,'Tol')
    MetaOpts.Optim.Tol = UQLabOpts.Optim.Tol;
end
if isfield(UQLabOpts.Optim,'Method')
    MetaOpts.Optim.Method = UQLabOpts.Optim.Method;
else
    MetaOpts.Optim.Method = 'HCMAES';
end
if isfield(UQLabOpts.Optim,'BFGS')
    MetaOpts.Optim.BFGS = UQLabOpts.Optim.BFGS;
end
if isfield(UQLabOpts.Optim,'GA')
    MetaOpts.Optim.GA = UQLabOpts.Optim.GA;
end
if isfield(UQLabOpts.Optim,'HGA')
    MetaOpts.Optim.HGA = UQLabOpts.Optim.HGA;
end
if isfield(UQLabOpts.Optim,'CMAES')
    MetaOpts.Optim.CMAES = UQLabOpts.Optim.CMAES;
end
if isfield(UQLabOpts.Optim,'HCMAES')
    MetaOpts.Optim.HCMAES = UQLabOpts.Optim.HCMAES;
end

% regression options
if Noise
    MetaOpts.Regression.SigmaNSQ = 'auto';
else
    MetaOpts.Regression.SigmaNSQ = 'none';
end

% scaling
MetaOpts.Scaling = false;

% train PCE-GPR model
M_gpr = uq_createModel(MetaOpts);


%% retrieve model parameters

rho = 10.^M_gpr.Kriging.theta; % kernel hyperparameter(s)
if Noise
    sigma_n = sqrt(M_gpr.Kriging.sigmaNSQ); % noise standard devation
else
    sigma_n = 0;
end
sigma = sqrt(M_gpr.Kriging.sigmaSQ); % kernel standard deviation
tau = sigma_n^2/(sigma^2+sigma_n^2); % noise ratio
Rtilde = M_gpr.Internal.Kriging.GP.R; % correlation matrix of training samples
L = chol(Rtilde+nugget*eye(size(Rtilde)),'lower'); % Cholesky factor of correlation matrix
alp = (1-tau)*(L'\(L\y_train)); % GPR model coefficients


%% assemble output structure

M_pce_gpr.x_train = x_train; % note: stored training samples are rescaled to standardized distribution
M_pce_gpr.y_train = y_train;
M_pce_gpr.distr = distr;
M_pce_gpr.basisfun = basis_functions;
M_pce_gpr.kernel = @(x,xp) kernel(x,xp,rho);
M_pce_gpr.L = L;
M_pce_gpr.rho = rho;
M_pce_gpr.tau = tau;
M_pce_gpr.sigma = sigma;
M_pce_gpr.sigma_n = sigma_n;
M_pce_gpr.alpha = alp;

% UQLab class model
M_pce_gpr.model = M_gpr;