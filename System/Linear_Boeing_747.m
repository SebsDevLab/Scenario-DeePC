%% File Name: Linear_Boeing_747.m
% Author: Sebastian Zieglmeier
% Date last updated: 10.06.2026
% Description: Linear (LTI) discrete-time state-space model of the Boeing 747
% longitudinal/lateral dynamics, used as a benchmark system for DeePC and
% Scenario-DeePC.
% Sources:
% [1] - Sebastian Zieglmeier, Nikolas Recke, Mathias Hudoba de Badyn,
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
function sys = Linear_Boeing_747()
    sys.name = "Linear_Boeing_747";
    %% Linear Boeing 747 state-space model
    A = [0.9997 0.0038 -0.0001 -0.0322;
        -0.0056 0.9648 0.7446 0.0001;
        0.0020 -0.0097 0.9543 -0.0000;
        0.0001 -0.0005 0.0978 1.0000];
    B = [0.0010 0.1000;
        -0.0615 0.0183;
        -0.1133 0.0586;
        -0.0057 0.0029];
    C = [1.0000 0 0 0;
        0 -1.0000 0 7.7400];
    D = [0];
    
    sys.T_samp = 0.1;

    %% LTI state-space model used for data collection
    sys.model.A = A;
    sys.model.B = B;    
    sys.model.C = C;    
    sys.model.D = D;


    %% Linearized state-space model used as the parametric model
    
    sys.param_model.A_M = A;
    sys.param_model.B_M = B;
    sys.param_model.C_M = C;
    sys.param_model.D_M = D;

    %% LTI state-space model used for simulation
    % sim_model_1: nominal scenario  -> data-collection system == simulation system
    % sim_model_2: robust scenario   -> data-collection system ~= simulation system

    sys.sim_model_1.A_sim = sys.model.A;
    sys.sim_model_1.B_sim = sys.model.B;
    sys.sim_model_1.C_sim = sys.model.C;
    sys.sim_model_1.D_sim = sys.model.D;


    %% System dimensions
    sys.nx = size(sys.model.A, 1);      % Number of states
    sys.nu = size(sys.model.B, 2);      % Number of inputs
    sys.ny = size(sys.model.C, 1);      % Number of outputs

    sys.nx_M = size(sys.param_model.A_M, 1); % Number of states of the parametric model
    sys.nu_M = size(sys.param_model.B_M, 2); % Number of inputs of the parametric model
    sys.ny_M = size(sys.param_model.C_M, 1); % Number of outputs of the parametric model
    
    %% Constraints
    sys.constraints.u_min = -20;
    sys.constraints.u_max = 20;
    sys.constraints.x_min = -inf;
    sys.constraints.x_max = inf;
    sys.constraints.y_min = [-25, -15];
    sys.constraints.y_max = [25, 15];
        
end