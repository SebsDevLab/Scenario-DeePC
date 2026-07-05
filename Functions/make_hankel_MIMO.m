%% File Name: make_hankel_MIMO.m
% Author: Sebastian Zieglmeier
% Date last updated: 10.06.2026
% Description:
% Constructs a block Hankel matrix for multi-input multi-output (MIMO) data
% sequences. This function reshapes and stacks the input data into a Hankel
% structure suitable for use in Data-Enabled Predictive Control (DeePC) and
% other data-driven control formulations.
%
% Usage:
%   H = make_hankel_MIMO(x, num_hankel_cols, L)
%
% Inputs:
%   x               : Measured input or output data sequence (m x K)
%   num_hankel_cols : Number of columns in the resulting Hankel matrix.
%   L               : Number of block rows (corresponding to the window length).
%
% Outputs:
%   H               : Constructed Hankel matrix of size (m*L x num_hankel_cols).
%
% Notes:
% 

function H = make_hankel_MIMO(x, num_hankel_cols, L)
    m = size(x, 1);
    H = zeros(m*L, num_hankel_cols);

    for i = 1:num_hankel_cols
        H(:, i) = reshape(x(:, i:i+L-1), [], 1);
    end
end