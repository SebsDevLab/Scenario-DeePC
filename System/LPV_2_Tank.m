%% File Name: LPV_2_Tank.m
% Author: Sebastian Zieglmeier
% Date last updated: 10.06.2026
% Description: Nonlinear cascaded two-tank system modeled as an LPV
% state-space system, used as a benchmark for DeePC and Scenario-DeePC.
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
%   (none)
%
% Outputs:
%   sys: struct holding the system matrices, dimensions, and constraints
%
% Notes:
%


function sys = LPV_2_Tank()
    sys.name = "LPV_2_Tank";
    %% Nonlinear two-tank LPV state-space model from [1]
    A_theta = @(theta_1, theta_2) [-0.904 * theta_1, 0; 
                                    0.904 * theta_1, -0.508 * theta_2];
    B = [0.258;
         0];
    C = [0,1];
    D = [0];
    
    sys.model.A = A_theta;
    sys.model.B = B;
    sys.model.C = C;
    sys.model.D = D;
    
    %% Compute the sampling time at a chosen linearization point
    x_1 = 1;    % Water level tank 1
    x_2 = 1;    % Water level tank 2
    theta_1 = 1/sqrt(x_1);
    theta_2 = 1/sqrt(x_2);
    A = A_theta(theta_1, theta_2);

    sys_c = ss(A, B, C, D);
    poles = eig(sys_c);     % Poles of the continuous-time system
    tau = 1./poles;
    tau_fast = min(abs(tau));   % Fastest time constant
    
    f_tau = 1/(2*pi*tau_fast); 
    f_samp = 10*f_tau;          % Sample at 10x the fastest dynamics

    T_samp = 1/f_samp;
    sys.T_samp = T_samp;

    %% Linearized, discretized state-space model used as the parametric model
    x_1 = 5;
    x_2 = 15;
    theta_1 = 1/sqrt(x_1);
    theta_2 = 1/sqrt(x_2);

    A = A_theta(theta_1, theta_2);
    sys2 = ss(A, B, C, D);
    sys_d = c2d(sys2, T_samp);

    sys.param_model.A_M = sys_d.A;
    sys.param_model.B_M = sys_d.B;
    sys.param_model.C_M = sys_d.C;
    sys.param_model.D_M = sys_d.D;

    %% LPV state-space model used for simulation
    % sim_model_1: nominal scenario  -> data-collection system == simulation system
    % sim_model_2: robust scenario   -> data-collection system ~= simulation system

    sys.sim_model_1.A_sim = sys.model.A;
    sys.sim_model_1.B_sim = sys.model.B;
    sys.sim_model_1.C_sim = sys.model.C;
    sys.sim_model_1.D_sim = sys.model.D;
    
    sys.sim_model_2.A_sim = @(theta_1,theta_2) 0.9 * A_theta(theta_1, theta_2);
    sys.sim_model_2.B_sim = sys.model.B;
    sys.sim_model_2.C_sim = sys.model.C;
    sys.sim_model_2.D_sim = sys.model.D;

    %% System dimensions
    sys.nx = size(sys.model.A(1,1), 1); % Number of states
    sys.nu = size(sys.model.B, 2);      % Number of inputs
    sys.ny = size(sys.model.C, 1);      % Number of outputs

    sys.nx_M = size(sys.param_model.A_M, 1); % Number of states of the parametric model
    sys.nu_M = size(sys.param_model.B_M, 2); % Number of inputs of the parametric model
    sys.ny_M = size(sys.param_model.C_M, 1); % Number of outputs of the parametric model

    
    
    %% Constraints
    sys.constraints.u_min = 0;
    sys.constraints.u_max = 22;
    sys.constraints.x_min = 0;
    sys.constraints.x_max = 100;
    sys.constraints.y_min = 0;
    sys.constraints.y_max = 25;

end