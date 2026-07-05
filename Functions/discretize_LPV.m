%% File Name: discretize_LPV.m
% Author: Sebastian Zieglmeier
% Date last updated: 10.06.2026
% Description: Discretizes the continuous-time LPV two-tank model at a given
% operating point so it can be simulated step by step.
% Sources:
% [1] - Teppa-Garran, et al., "Liquid level tracking for a coupled tank system
%       using quasi-LPV control.", Ingenius 33 (2025): 15-26.
% [2] - Sebastian Zieglmeier, Nikolas Recke, Mathias Hudoba de Badyn,
%       "Scenario-based Data-Enabled Predictive Control: Robustification via
%       the Scenario Approach".
%       Code: https://github.com/SebsDevLab/Scenario-DeePC.git
%       (TODO: add arXiv/DOI once available)
%
% Inputs:
%   A_theta, B_c, C_c, D_c: matrices of the continuous-time LPV system
%                           (A_theta is a function handle of the scheduling parameters)
%   x: current state, used to evaluate the scheduling parameters
%   T_samp: sampling time for discretization
% Outputs:
%   A, B, C, D: discrete-time system matrices
%
% Notes:
%


function [A, B, C, D] = discretize_LPV(A_theta, B_c, C_c, D_c, x, T_samp)
    % Guard against division by zero in the scheduling parameters
    x_1 = max(x(1,1), 1);
    x_2 = max(x(1,2), 1);
    theta_1 = 1/sqrt(x_1);
    theta_2 = 1/sqrt(x_2);
    % Evaluate the LPV system at this operating point and discretize
    A = A_theta(theta_1, theta_2);
    sys2 = ss(A, B_c, C_c, D_c);
    sys_d = c2d(sys2, T_samp);
    A = sys_d.A;
    B = sys_d.B; 
    C = sys_d.C; 
    D = sys_d.D;
end
    