function indices = PCEindex(p, d, u)

% indices = PCEindex(p, d, u) generates a matrix with the d-dimensional set
% of PCE multi-indices of total degree bounded by p. The multi-indices are
% ordered by increasing total degree. The optional parameter 0 < u < 1
% defines a hyperbolic truncation that bounds the u-norm of the
% multi-indices by p. Providing no u is equivalent to setting u = 1.
%
% Inputs:
%   p - maximum total degree (non-negative integer)
%   d - number of dimensions (positive integer)
%   u - norm index for hyperbolic truncation (0 < u < 1)
%
% Output:
%   indices - matrix where each row represents a multi-index
%
% Acknowledgment:
% This code was developed with assistance from OpenAI's ChatGPT (GPT-4),
% and then customized.
%
% Author: Paolo Manfredi (with the help of ChatGPT - GPT-4)
% Affiliation: Politecnico di Torino
% Date: July 2024

    % Calculate the total number of multi-indices
    numIndices = nchoosek(p + d, d);
    
    % Preallocate the result matrix
    indices = zeros(numIndices, d);
    
    % Index to keep track of the position in the result matrix
    currentIndex = 1;
    
    % Loop through each total degree from 0 to p
    for total_degree = 0:p
        % Generate all combinations of degrees for the current total degree
        combs = generateCombinations(total_degree, d);
        numCombs = size(combs, 1);
        
        % Fill in the preallocated matrix with the combinations
        indices(currentIndex:currentIndex + numCombs - 1, :) = combs;
        currentIndex = currentIndex + numCombs;
    end

    % Hyperbolic truncation
    if nargin>2
        % Remove elements with hyperolic degree larger than p
        indices(vecnorm(indices,u,2)>p,:) = [];
    end
end

function combs = generateCombinations(total_degree, d)
    % This function generates all combinations of d non-negative integers
    % that sum to total_degree.
    
    if d == 1
        combs = total_degree;
        return;
    end
    
    % Calculate the number of combinations for preallocation
    numCombs = nchoosek(total_degree + d - 1, d - 1);
    
    % Preallocate the combinations matrix
    combs = zeros(numCombs, d);
    
    % Index to keep track of the position in the combinations matrix
    currentIndex = 1;
    
    for k = 0:total_degree
        sub_combs = generateCombinations(total_degree - k, d - 1);
        numSubCombs = size(sub_combs, 1);
        
        % Fill in the preallocated matrix with the combinations
        combs(currentIndex:currentIndex + numSubCombs - 1, :) = [k * ones(numSubCombs, 1), sub_combs];
        currentIndex = currentIndex + numSubCombs;
    end
end
