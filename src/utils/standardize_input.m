function x_standard = standardize_input(x,distr)

% x_standard = standardize_input(x,distr) rescale input values x to
% the standardized version of the distribution specified by distr, i.e.,
% N(0,1) for the Gaussian distribution and U(-1,1) for the uniform
% distribution
%
%   x: matrix of input samples of size n x d, where n is the number of
%   samples and d is the number of dimensions
%
%   distr: structure array describing the probability distributions of
%   input parameters (see train_pce_gpr)
%
% Author: Paolo Manfredi
% Affiliation: Politecnico di Torino
% Date: February 2025

[N,d] = size(x);
x_standard = zeros(N,d);

for jj = 1:d
    switch lower(distr(jj).Type)
        case 'norm'
            m = distr(jj).Parameters(1);
            s = distr(jj).Parameters(2);
            x_standard(:,jj) = (x(:,jj) - m)/s;

        case 'unif'
            a = distr(jj).Parameters(1);
            b = distr(jj).Parameters(2);
            x_standard(:,jj) = (x(:,jj) - (a+b)/2)/((b-a)/2);
    end
end