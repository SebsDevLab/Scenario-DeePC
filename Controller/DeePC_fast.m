%% File Name: DeePC_fast.m
% Author: Sebastian Zieglmeier
% Date last updated: 10.06.2026
% Description: Accelerated standard DeePC controller. The optimization problem
% is precompiled once as a YALMIP parametric solver and only evaluated at
% runtime. Used as the nominal baseline controller for the Scenario-DeePC study.
%
% Sources:
% [1] - Sebastian Zieglmeier, Nikolas Recke, Mathias Hudoba de Badyn,
%       "Scenario-based Data-Enabled Predictive Control: Robustification via
%       the Scenario Approach".
%       Code: https://github.com/SebsDevLab/Scenario-DeePC.git
%       (TODO: add arXiv/DOI once available)
% [2] - Sebastian Zieglmeier, et al., "Data-Enabled Predictive Control and
%       Guidance for Autonomous Underwater Vehicles",
%       https://doi.org/10.48550/arXiv.2510.25309 (original fast-DeePC implementation)
%
% Properties:
% T_d: Number of data collection points
% T_ini: Number of initialization points
% T_fut: Prediction horizon
% R, Q: Cost matrices for input and output
% lambda_ini, lambda_g: Regularization hyperparameters
% H_u, H_y: Hankel/Page matrices for input and output data
% u_min, u_max: Input constraints
% y_min, y_max: Output constraints
%
% Internal variables:
% g, sigma, u, y: Decision variables of the DeePC optimization problem
% u_past_par, y_past_par, y_ref_par: Parametric placeholders for past data and reference
% paramSolver: Precompiled parametric solver object (YALMIP 'optimizer')
% options: Solver options (e.g., MOSEK, verbosity level)
%
% Inputs (step function):
% u_past, y_past: Past T_ini input and output data
% y_ref: Reference trajectory over the prediction horizon
%
% Outputs (step function):
% u_fut: Computed control input sequence
% y_value, g_value, sigma_value: Optimization results for analysis
% yalmip_error_flag: Error flag (0 if successful, 1 if solver issue)
%
% 
% Motivation:
% The classic DeePC formulation repeatedly constructs and solves a large-scale 
% optimization problem at every control step, which causes high computational 
% load. This implementation uses YALMIP's 'optimizer' object to precompile the 
% optimization problem once, treating the past data and reference trajectory 
% as parameters. At runtime, only the parametric evaluation is performed, 
% which drastically reduces execution time while preserving identical 
% control behavior.
% Further Notes:
% - This implementation achieves the same control law as the standard DeePC
%   but runs significantly faster, making it suitable for real-time and
%   embedded applications.
% - The problem structure (constraints and cost) remains identical to the
%   original DeePC formulation.
% - Compilation is done only once via initializeOptimization(); subsequent
%   calls only perform fast parametric evaluation.



classdef DeePC_fast < handle
    properties
        % Problem dims & weights
        T_d; T_ini; T_fut; T;
        R; Q; lambda_ini; lambda_g;
        nu; ny;
        H_u; H_y;

        % bounds
        u_min; u_max; y_min; y_max;

        % decision vars
        g; sigma; u; y;

        % parameter placeholders (sdpvar shapes)
        u_past_par; y_past_par; y_ref_par;

        % compiled parametric solver
        paramSolver;    % optimizer object

        % solver options flag
        initialized = false;
        options;
    end

    methods
        function obj = DeePC_fast(T_d, T_ini, T_fut, R, Q, lambda_ini, lambda_g, H_u, H_y, sys)
            obj.T_d = T_d; obj.T_ini = T_ini; obj.T_fut = T_fut;
            obj.T = T_ini + T_fut;
            obj.R = R; obj.Q = Q;
            obj.lambda_ini = lambda_ini; obj.lambda_g = lambda_g;
            obj.H_u = H_u; obj.H_y = H_y;
            obj.nu = sys.nu; obj.ny = sys.ny;
            obj.u_min = sys.constraints.u_min .* ones(T_fut, sys.nu);
            obj.u_max = sys.constraints.u_max .* ones(T_fut, sys.nu);
            obj.y_min = sys.constraints.y_min .* ones(T_fut, sys.ny);
            obj.y_max = sys.constraints.y_max .* ones(T_fut, sys.ny);

            % Default solver options
            obj.options = sdpsettings('verbose', 0, 'solver', 'mosek', 'debug', 0);
        end

        function initializeOptimization(obj)
            if obj.initialized
                return
            end

            % --- Decision variables ---
            ncols_g = obj.T_d - obj.T + 1;
            if ncols_g <= 0
                error('T_d - T + 1 must be > 0.');
            end
            obj.g     = sdpvar(ncols_g, 1);
            obj.sigma = sdpvar(obj.T_ini, obj.ny);     % T_ini x ny
            obj.u     = sdpvar(obj.T_fut, obj.nu);     % T_fut x nu
            obj.y     = sdpvar(obj.T_fut, obj.ny);     % T_fut x ny

            % --- Parameter placeholders (runtime inputs to the optimizer) ---
            obj.u_past_par = sdpvar(obj.T_ini, obj.nu, 'full');   % T_ini x nu
            obj.y_past_par = sdpvar(obj.T_ini, obj.ny, 'full');   % T_ini x ny
            obj.y_ref_par  = sdpvar(obj.T_fut, obj.ny, 'full');   % T_fut x ny

            % --- Hankel sub-blocks (constant) ---
            U_p = obj.H_u(1:obj.nu * obj.T_ini, :);
            U_f = obj.H_u(obj.nu * obj.T_ini + 1:end, :);
            Y_p = obj.H_y(1:obj.ny * obj.T_ini, :);
            Y_f = obj.H_y(obj.ny * obj.T_ini + 1:end, :);

            % --- Cost: initial-condition slack + g-regularization ---
            cost = reshape(obj.sigma, [], 1)' * obj.lambda_ini * reshape(obj.sigma, [], 1) ...
                 + obj.lambda_g * (obj.g' * obj.g);
            
            % Tracking and input-effort cost over the prediction horizon
            for i = 1:obj.T_fut
                cost = cost + (obj.y(i,:) - obj.y_ref_par(i,:)) * obj.Q * (obj.y(i,:) - obj.y_ref_par(i,:))' ...
                             + obj.u(i,:) * obj.R * obj.u(i,:)';
            end

            % --- Constraints ---
            cons = [ ...
                U_p * obj.g == reshape(obj.u_past_par', [], 1);
                U_f * obj.g == reshape(obj.u', [], 1);
                Y_p * obj.g == reshape((obj.y_past_par + obj.sigma)', [], 1);
                Y_f * obj.g == reshape(obj.y', [], 1);
                obj.y(end,:) == obj.y_ref_par(end,:);
                obj.u(:) >= obj.u_min(:);
                obj.u(:) <= obj.u_max(:);
                obj.y(:) >= obj.y_min(:);
                obj.y(:) <= obj.y_max(:);
            ];

            % --- Compile the parametric solver (parameters: past I/O and reference) ---
            obj.paramSolver = optimizer(cons, cost, obj.options, ...
                                       {obj.u_past_par, obj.y_past_par, obj.y_ref_par}, ...
                                       {obj.u, obj.y, obj.g, obj.sigma});
            obj.initialized = true;
        end

        function [u_fut, y_value, g_value, sigma_value, yalmip_error_flag] = step(obj, u_past, y_past, y_ref, ~, ~)
            % The two ignored arguments (H_u, H_y) are kept for a common step()
            % interface with the other controllers; the Hankel data is already
            % baked into the compiled solver here.
            % Build/compile once if needed
            if ~obj.initialized
                obj.initializeOptimization();
            end

            % --- Call the compiled parametric solver ---
            [sol, diagnostics] = obj.paramSolver({u_past, y_past, y_ref});

            % diagnostics may be numeric (0) or a struct with field 'problem' depending on YALMIP build.
            ok = false;
            if isempty(diagnostics)
                ok = true;
            elseif isnumeric(diagnostics)
                ok = (diagnostics == 0);
            elseif isstruct(diagnostics) && isfield(diagnostics, 'problem')
                ok = (diagnostics.problem == 0);
            end

            if ok
                yalmip_error_flag = 0;
            else
                yalmip_error_flag = 1;
            end
            u_fut = sol{1};
            y_value = sol{2};
            g_value = sol{3};
            sigma_value = sol{4};
        end

        
    end
end
