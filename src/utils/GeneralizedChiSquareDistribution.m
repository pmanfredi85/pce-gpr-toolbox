classdef GeneralizedChiSquareDistribution

%GeneralizedChiSquareDistribution model
%   This is a GeneralizedChiSquareDistribution object. The class can be
%   used to compute PDF, CDF, inverse CDF, and to generate random samples
%   from the distribution
%
%   The class builds upon the Generalized chi-square distribution toolbox
%   by Abhranil Das - version 2.3.0
%
%   https://it.mathworks.com/matlabcentral/fileexchange/85028-generalized-chi-square-distribution
%
%   A. Das, "New methods for computing the generalized chi-square
%   distribution." arXiv preprint arXiv:2404.05062 (2024).
%
%   GeneralizedChiSquareDistribution properties:
%       w           - weights (row vector)
%       k           - degrees of freedom (row vector)
%       lambda      - non-centrality parameters (row vector)
%       m           - mean of normal term
%       s           - standard deviation of normal term
%
%   GeneralizedChiSquareDistribution methods:
%       pdf         - computes PDF
%       cdf         - computes CDF
%       icdf        - computes inverse CDF
%       ci          - computes confidence interval
%       random      - generates random samples
%       mean        - returns distribution mean
%       var         - returns distribution variance
%       std         - returns distribution standard deviation
%
% Author: Paolo Manfredi
% Affiliation: Politecnico di Torino
% Date: April 2024

    properties
        % W - weights (row vector)
        w;

        % K - degrees of freedom (row vector)
        k;

        % LAMBDA - non-centrality parameters (row vector)
        lambda;

        % M - mean of normal term
        m;

        % S - standard deviation of normal term
        s;

    end

    methods
        % Constructor
        function obj = GeneralizedChiSquareDistribution(w,k,lambda,m,s)
            obj.w = w(:);
            obj.k = k(:);
            obj.lambda = lambda(:);
            obj.m = m;
            obj.s = s;
        end

        % PDF method - Computes PDF at points specified by x. See gx2pdf
        % for optional arguments
        function f = pdf(this,x,varargin)
            f = gx2pdf(x,this.w',this.k',this.lambda',this.m,this.s,varargin{:});
        end

        % CDF method - Computes CDF at points specified by x. See gx2cdf
        % for optional arguments
        function p = cdf(this,x,varargin)
            p = gx2cdf(x,this.w',this.k',this.lambda',this.m,this.s,varargin{:});
        end

        % ICDF method - Computes inverse CDF at points specified by x. See
        % gx2inv for optional arguments
        function x = icdf(this,p,varargin)
            x = gx2inv(p,this.w',this.k',this.lambda',this.m,this.s,varargin{:});
        end

        % CI method - Computes confidence interval with significance level
        % alpha
        function c = ci(this,alpha)
            c = [this.icdf(alpha/2), this.icdf(1-alpha/2)];
        end

        % RANDOM method - Generates random samples of specified size. See
        % gx2rnd for details
        function r = random(this,varargin)
            r = gx2rnd(this.w',this.k',this.lambda',this.m,this.s,varargin{:});
        end

        % MEAN method - Returns distribution mean. See gx2stat for details
        function mu = mean(this)
            mu = gx2stat(this.w',this.k',this.lambda',this.m,this.s);
        end

        % VAR method - Returns distribution variance. See gx2stat for
        % details
        function v = var(this)
            [~,v] = gx2stat(this.w',this.k',this.lambda',this.m,this.s);
        end

        % STD method - Returns distribution standard deviation. See gx2stat
        % for details
        function st = std(this)
            [~,v] = gx2stat(this.w',this.k',this.lambda',this.m,this.s);
            st = sqrt(v);
        end

    end

end