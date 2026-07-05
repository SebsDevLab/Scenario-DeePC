%% File Name: simulate_Scenario_DeePC_Linear_Boeing_Online.m
% Author: Sebastian Zieglmeier
% Date last updated: 10.06.2026
% Description: Closed-loop simulation of the linear Boeing 747 model controlled by
% DeePC and Scenario-DeePC, with an ONLINE (adaptive) scenario buffer that is
% updated during operation so the uncertainty description tracks time-varying
% noise and disturbances. Standard DeePC is run first to seed the buffer, then
% adaptive Scenario-DeePC is run for comparison.
% Sources:
% [1] - Sebastian Zieglmeier, Nikolas Recke, Mathias Hudoba de Badyn,
%       "Scenario-based Data-Enabled Predictive Control: Robustification via
%       the Scenario Approach".
%       Code: https://github.com/SebsDevLab/Scenario-DeePC.git
%       (TODO: add arXiv/DOI once available)
%
% Notes:
% 


close all;
clc;
clear all;


%% Get system
sys_name = "Linear_Boeing_747";
sys = eval(sys_name);

%% Select active noise mode
active_noise = 3;   
% 1 - no noise
% 2 - measurement + process noise
% 3 - measurement noise only
% 4 - process noise only

%% Get controller
control_name = "DeePC_fast";   % options: DeePC_fast, MPC

%% Initialize variables and control hyperparameters
T_sim = 400;       % Simulation length [steps]
T_fut = 20;         % Prediction horizon N
T_d = 1000;         % Number of data points used to build the Hankel matrices
T_ini = 20;         % Length of the initial (past) window T_ini

% DeePC cost-function weights and regularization parameters
lambda_ini = 1e7;   % Penalty on the initial-condition slack sigma_y
lambda_g = 1e5;     % Penalty on the regularizer ||g||^2
r = 1e-3;           % Input weight R
q = 1e2;            % Output (tracking) weight Q

%% Load system
% Linear Boeing 747 state-space model used to generate the data
A = sys.model.A;
B = sys.model.B;
C = sys.model.C;
D = sys.model.D;

A_sim = sys.sim_model_1.A_sim;
B_sim = sys.sim_model_1.B_sim;
C_sim = sys.sim_model_1.C_sim;
D_sim = sys.sim_model_1.D_sim;


% System constants
nx = sys.nx; % Number of states
nu = sys.nu; % Number of inputs
ny = sys.ny; % Number of outputs

constraints = sys.constraints;

%% Collect data for the data-driven controller (DeePC)

u_data = zeros(T_d, nu);
u_data_sys = u_data;   

y_data = zeros(T_d, ny);
y_data_sys = y_data;
y_meas = y_data;
v_meas = y_data;
v_proc = zeros(T_d, nx);

x_data_sys = zeros(T_d+1, nx);   

% NOTE: v_meas / v_proc here are the raw measurement and process noise
% sequences used during data collection, NOT the papers prediction-error
% scenarios w^{(j)} (those are built later from the prediction error).
sigma_meas = 0.01 * sys.constraints.y_max;  % measurement-noise std, scaled to max output
sigma_proc = 0.1*sigma_meas;                % process-noise std

rng(1,'twister');
%rand('seed', 8); % seeding for reproducibility
u_prbs = 1.5*idinput([T_d,2], 'RBS');
for i = 1:T_d
    randn(1,nu);                          % advance RNG to reproduce published noise stream
    u_data_sys(i,:) = u_prbs(i,:);
    [u_data_sys(i, :), w, warn] = system_boundaries(u_data_sys(i, :), sys.constraints, "u");
    if w == 1
        disp(warn);
    end

    if active_noise == 2
        v_meas(i,:) = sigma_meas .* randn(1, ny);
        v_proc(i,:) = sigma_proc .* randn(1, nx);
    elseif active_noise == 3
        v_meas(i,:) = sigma_meas .* randn(1, ny);
        v_proc(i,:) = zeros(1, nx);
    elseif active_noise == 4
        v_meas(i,:) = zeros(1, ny);
        v_proc(i,:) = sigma_proc .* randn(1, nx);
    else
        v_meas(i,:) = zeros(1, ny);
        v_proc(i,:) = zeros(1, nx);
    end
    x_data_sys(i+1, :) = (A * x_data_sys(i, :)' + B * u_data_sys(i,:)')' + v_proc(i,:); 
    y_data_sys(i,:) = (C * x_data_sys(i, :)')'; 
    y_meas(i,:) = y_data_sys(i,:) + v_meas(i,:);
    
    % Enforce state and output constraints on the collected data
    [x_data_sys(i+1, :), w, warn] = system_boundaries(x_data_sys(i+1, :), sys.constraints, "x");
    if w == 1
        disp(warn);
    end
    [y_data_sys(i, :), w, warn] = system_boundaries_multi(y_data_sys(i, :), sys.constraints, "y");
    if w >= 1
        disp(warn);
    end
end
u_data = u_data_sys;
y_data = y_meas;

%% Diagnostics: collected data and measurement-noise distribution
figure
plot(y_data_sys);
hold on;
plot(y_data);
legend on;

w = v_meas(:);
mu_hat = mean(w);
sigma_hat = std(w);
figure
histogram(w, 'Normalization', 'pdf');
hold on;
xline(mu_hat, "r--", 'LineWidth', 2);
% x-grid for the fitted Gaussian
x_vals = linspace(min(w), max(w), 200);

% Gaussian PDF
pdf_vals = (1/(sigma_hat*sqrt(2*pi))) * ...
           exp(-0.5*((x_vals - mu_hat)/sigma_hat).^2);

% Plot fitted Gaussian
plot(x_vals, pdf_vals, 'LineWidth', 2);





%% Build Hankel matrices
L = T_ini + T_fut; % Window depth L = T_ini + N (sets the Hankel block height)
num_hankel_cols = T_d - L +1;
H_u = make_hankel_MIMO(u_data', num_hankel_cols, L);
H_y = make_hankel_MIMO(y_data', num_hankel_cols, L);




%% Initialize controller and simulation variables

if control_name == "MPC"
    control = MPC(T_d, T_ini, T_fut, r, q, lambda_ini, lambda_g, H_u, H_y, sys);
elseif control_name == "DeePC_fast"
    control = DeePC_fast(T_d, T_ini, T_fut, r, q, lambda_ini, lambda_g, H_u, H_y, sys);
else
    disp("Warning: No controller selected. Check spelling of control_name.");
end

% Preallocate arrays to store simulation results
u_sim = zeros(T_sim, nu);
x_sim = zeros(T_sim + 1, nx);
y_sim = zeros(T_sim, ny);
y_sim_pred = zeros(T_sim, ny);
y_sim_meas = zeros(T_sim, ny);
w_sim_meas3 = zeros(T_sim, ny);
w_sim_proc3 = zeros(T_sim, nx);

y_D = zeros(T_sim, ny);

%% Initial condition
x_sim(1, :) = zeros(1, nx);
u_past_sim = zeros(T_ini, nu);
y_past_sim = zeros(T_ini, ny);

%% Reference trajectory

ref_name = "Smooth_Step";
ref(1,:) = 2*get_ref(ref_name, T_sim, 1, T_fut, T_ini);
ref(2,:) = -1*get_ref(ref_name, T_sim, 1, T_fut, T_ini);
% Available reference profiles: Smooth_Step_to_Sinus, Smooth_Step, Sinus, Step
% Reference values can be changed in get_ref.m
over_step = 0;
%% Simulation loop
for i = 1:T_sim 
    y_reference = ref(:,i:i+T_fut-1);
    [u_fut, y_fut, g_fut, sigma_fut, opt_error_flag] = control.step(u_past_sim, y_past_sim, y_reference', H_u, H_y); 
    if opt_error_flag == 1
        disp(["Warning: Opt_error, Timestep: " + string(i)]);
    end
    u_sim(i,:) = u_fut(1,:); % Apply only the first input (receding horizon)
    % Enforce input constraints
    [u_sim(i, :), w, warn] = system_boundaries(u_sim(i, :), sys.constraints, "u");
    if w == 1
        disp(warn);
    end

    if active_noise == 2
        w_sim_meas3(i,:) = sigma_meas .* randn(1, ny);
        w_sim_proc3(i,:) = sigma_proc .* randn(1, nx);
    elseif active_noise == 3
        w_sim_meas3(i,:) = sigma_meas .* randn(1, ny);
        w_sim_proc3(i,:) = zeros(1, nx);
    elseif active_noise == 4
        w_sim_meas3(i,:) = zeros(1, ny);
        w_sim_proc3(i,:) = sigma_proc .* randn(1, nx);
    else
        w_sim_meas3(i,:) = zeros(1, ny);
        w_sim_proc3(i,:) = zeros(1, nx);
    end
    % Simulate the real system
    x_sim(i+1, :) = (A_sim * x_sim(i, :)' + B_sim * u_sim(i,:)')' + w_sim_proc3(i,:);
    y_sim(i, :) = (C_sim * x_sim(i, :)')';
    y_sim_pred(i, :) = y_fut(1,:);
    y_sim_meas(i, :) = y_sim(i,:) + w_sim_meas3(i,:);

    
    % Enforce state constraints
    [x_sim(i+1, :), w, warn] = system_boundaries(x_sim(i+1, :), sys.constraints, "x"); 
    if w == 1
        disp(warn);
    end
    [~, w, warn] = system_boundaries_multi(y_sim_meas(i, :), sys.constraints, "y"); 
    if w >= 1
        over_step = over_step + w;   % count output-constraint violations
    end

    % Store the measured output for the past window
    y_D(i,:) = y_sim_meas(i,:);
    % Update the past input/output windows for the next step
    u_past_sim = [u_past_sim(2:end,:); u_sim(i,:)];
    y_past_sim = [y_past_sim(2:end,:); y_D(i,:)];
end
disp("Constraint violations without Scenario approach: ");
disp(over_step);
e_2 = y_sim_meas - y_sim_pred;        % measured prediction error w^{(j)} (used for scenarios)

% Preallocate arrays to store simulation results
u_sim = zeros(T_sim, nu);
x_sim = zeros(T_sim + 1, nx);
y_sim = zeros(T_sim, ny);
y_sim_pred = zeros(T_sim, ny);
y_sim_meas = zeros(T_sim, ny);
w_sim_meas = zeros(T_sim, ny);
w_sim_proc = zeros(T_sim, nx);

y_D = zeros(T_sim, ny);

%% Initial condition
x_sim(1, :) = zeros(1, nx);
u_past_sim = zeros(T_ini, nu);
y_past_sim = zeros(T_ini, ny);

%% Reference trajectory

ref_name = "Smooth_Step";
ref(1,:) = 1.5*get_ref(ref_name, T_sim, 1, T_fut, T_ini);
ref(1,201:T_sim+T_fut) = 25 * ones(1,T_sim+T_fut-200);
ref(1, 5:end) = 25;
ref(2,:) = -1*get_ref(ref_name, T_sim, 1, T_fut, T_ini);
ref(2,201:T_sim+T_fut) = -10 * ones(1,T_sim+T_fut-200);
ref(2, 5:end) = -10;
% Available reference profiles: Smooth_Step_to_Sinus, Smooth_Step, Sinus, Step
% Reference values can be changed in get_ref.m
over_step = 0;
%% Simulation loop
for i = 1:T_sim 
    y_reference = ref(:,i:i+T_fut-1);
    [u_fut, y_fut, g_fut, sigma_fut, opt_error_flag] = control.step(u_past_sim, y_past_sim, y_reference', H_u, H_y); 
    if opt_error_flag == 1
        disp(["Warning: Opt_error, Timestep: " + string(i)]);
    end
    u_sim(i,:) = u_fut(1,:); % Apply only the first input (receding horizon)
    % Enforce input constraints
    [u_sim(i, :), w, warn] = system_boundaries(u_sim(i, :), sys.constraints, "u");
    if w == 1
        disp(warn);
    end

    if active_noise == 2
        w_sim_meas(i,:) = sigma_meas .* randn(1, ny);
        w_sim_proc(i,:) = sigma_proc .* randn(1, nx);
    elseif active_noise == 3
        % Time-varying scenario (motivates the online/adaptive buffer):
        % a slowly growing output disturbance switches on at step 300 ...
        if i >= 300
            disturbance = disturbance + [0.05,0];
            disturbance = min(disturbance, [1,0]);
        else
            disturbance = [0,0];
        end
        % ... while the measurement-noise level changes across three phases
        if i>= 200
            sigma_meas = sigma_meas + 0.001*sys.constraints.y_max;
            sigma_meas = min(sigma_meas, 0.01*sys.constraints.y_max);
        elseif i >= 100 
            sigma_meas = sigma_meas - 0.001*sys.constraints.y_max;
            sigma_meas = max(sigma_meas, 0.005*sys.constraints.y_max);
        else
            sigma_meas = 0.02 * sys.constraints.y_max;
        end
        w_sim_meas(i,:) = sigma_meas .* randn(1, ny) + disturbance;
        w_sim_proc(i,:) = zeros(1, nx);
    elseif active_noise == 4
        w_sim_meas(i,:) = zeros(1, ny);
        w_sim_proc(i,:) = sigma_proc .* randn(1, nx);
    else
        w_sim_meas(i,:) = zeros(1, ny);
        w_sim_proc(i,:) = zeros(1, nx);
    end
    % Simulate the real system
    x_sim(i+1, :) = (A_sim * x_sim(i, :)' + B_sim * u_sim(i,:)')' + w_sim_proc(i,:);
    y_sim(i, :) = (C_sim * x_sim(i, :)')';
    y_sim_pred(i, :) = y_fut(1,:);
    y_sim_meas(i, :) = y_sim(i,:) + w_sim_meas(i,:);

    
    % Enforce state constraints
    [x_sim(i+1, :), w, warn] = system_boundaries(x_sim(i+1, :), sys.constraints, "x"); 
    if w == 1
        disp(warn);
    end
    [~, w, warn] = system_boundaries_multi(y_sim_meas(i, :), sys.constraints, "y"); 
    if w >= 1
        over_step = over_step + w;   % count output-constraint violations
    end

    % Store the measured output for the past window
    y_D(i,:) = y_sim_meas(i,:);
    % Update the past input/output windows for the next step
    u_past_sim = [u_past_sim(2:end,:); u_sim(i,:)];
    y_past_sim = [y_past_sim(2:end,:); y_D(i,:)];
end
disp("Constraint violations without Scenario approach: ");
disp(over_step);

%% Numerical evaluation
y_ref = ref(:, 1:T_sim)';
rel_error= zeros(size(y_ref));
% Relative tracking error [%]
for i = 1:length(y_ref)
    for j = 1:ny
        if y_ref(i,j) ~= 0
            rel_error(i,j) = abs(y_sim(i,j) - y_ref(i,j)) / y_ref(i,j) * 100;
        else
            % Guard against division by zero when the reference is zero
            rel_error(i,j) = 0;
        end
    end
end
Avg_rel_error = mean(rel_error);
for j=1:ny
    RMSE(j) = 1/1 * rmse(ref(j,1:T_sim), y_sim_meas(1:T_sim,j)');
    % RMSE2(j) = 1/1 * rmse(ref(j,201:T_sim), y_sim_meas(201:T_sim,j)');
end
disp(RMSE)
% disp(RMSE2)
%% Graphic evaluation
figure
plot(0:T_sim-1, y_sim, linewidth=1.5)
hold on
plot(0:T_sim-1, ref(:,1:T_sim)', 'g--', linewidth=1)
xlabel('Timestep [-]')
ylabel('Height [cm]')
title("Output", 'FontSize', 14)
legend("DeePC", "reference");
set(gcf, 'Color', 'w');
grid on
exportgraphics(gcf, 'DeePC_general.png', ...
    'Resolution', 600);


figure
subplot(3, 1, 1)
plot(0:T_sim-1, y_sim_meas, linewidth=1.5)
hold on
plot(0:T_sim-1, ref(:,1:T_sim)', 'g--', linewidth=1)
xlabel('Discrete timestep')
title("Output", 'FontSize', 14)
legend("y", "y_{ref}");
grid on

subplot(3, 1, 2)
plot(0:T_sim-1, u_sim, 'g')
xlabel('Discrete timestep')
title("Control input", 'FontSize', 14)
grid on

subplot(3, 1, 3)
plot(0:T_sim-1, rel_error, 'r')
xlabel('Discrete timestep')
title("Relative Error", 'FontSize', 14)
grid on





% -------------------------------------------------------------------------
%% Scenario-DeePC run
% -------------------------------------------------------------------------

%% Initialize controller and simulation variables
N_scen = 25;
control = Scenario_DeePC_multi_Sc_cost(T_d, T_ini, T_fut, r, q, lambda_ini, lambda_g, 1e6, H_u, H_y, sys, N_scen);


% Preallocate arrays to store simulation results
u_sim2 = zeros(T_sim, nu);
x_sim2 = zeros(T_sim + 1, nx);
y_sim2 = zeros(T_sim, ny);
y_sim_pred2 = zeros(T_sim, ny);
y_sim_meas2 = zeros(T_sim, ny);
w_sim_meas2 = w_sim_meas;       % reuse identical noise realization for a fair comparison
w_sim_proc2 = w_sim_proc;       % reuse identical noise realization for a fair comparison
error_buffer = zeros(50, ny);   % scenario buffer E_buffer (filled online during the run)
y_D = zeros(T_sim, ny);

%% Initial condition
x_sim2(1, :) = zeros(1, nx);
u_past_sim = zeros(T_ini, nu);
y_past_sim = zeros(T_ini, ny);

%% Reference trajectory

% Reference trajectory is reused from the DeePC run above
disp("Scenario approach start");
over_step = 0;
%% Simulation loop (online Scenario-DeePC)
for i = 1:T_sim 
    y_reference = ref(:,i:i+T_fut-1);
    % Sample N_scen scenarios over the horizon from the current buffer E_buffer
    idx = randi(size(error_buffer, 1), T_fut * N_scen, 1);
    for j = 1:ny
        w_flat(:,j) = error_buffer(idx, j);
        w_scen(:,j,:) = reshape(w_flat(:,j), T_fut, 1, N_scen);
    end
    [u_fut, y_fut, g_fut, sigma_fut, opt_error_flag] = control.step(u_past_sim, y_past_sim, y_reference', H_u, H_y, w_scen); 
    if opt_error_flag == 1
        disp(["Warning: Opt_error, Timestep: " + string(i)]);
    end
    u_sim2(i,:) = u_fut(1,:); % Apply only the first input (receding horizon)
    % Enforce input constraints
    [u_sim2(i, :), w, warn] = system_boundaries(u_sim2(i, :), sys.constraints, "u");
    if w == 1
        disp(warn);
    end

    % Noise is intentionally NOT regenerated here: w_sim_meas2/w_sim_proc2 were
    % set to the same realizations used in the DeePC run above, so both
    % controllers see identical noise and the comparison is fair.

    % Simulate the real system
    x_sim2(i+1, :) = (A_sim * x_sim2(i, :)' + B_sim * u_sim2(i,:)')' + w_sim_proc2(i,:);
    y_sim2(i, :) = (C_sim * x_sim2(i, :)')';
    y_sim_pred2(i, :) = y_fut(1,:);
    y_sim_meas2(i, :) = y_sim2(i,:) + w_sim_meas2(i,:);
    
    

    % Enforce state constraints
    [x_sim2(i+1, :), w, warn] = system_boundaries(x_sim2(i+1, :), sys.constraints, "x"); 
    if w == 1
        disp(warn);
    end
    [~, w, warn] = system_boundaries_multi(y_sim_meas2(i, :), sys.constraints, "y"); 
    if w >= 1
        over_step = over_step + w;   % count output-constraint violations
    end

    % Store the measured output for the past window
    y_D(i,:) = y_sim_meas2(i,:);
    % Update the past input/output windows for the next step
    u_past_sim = [u_past_sim(2:end,:); u_sim2(i,:)];
    y_past_sim = [y_past_sim(2:end,:); y_D(i,:)];
    y_err(i,:) = y_sim_meas2(i,:) - y_sim_pred2(i,:);   % one-step prediction error w^{(j)}
    % Online (adaptive) update: push the newest error into the buffer every step
    % so the scenario set tracks the time-varying noise and disturbance.
    error_buffer = [error_buffer(2:end, :); y_err(i,:)];
end
disp("Constraint violations with Scenario approach: ");
disp(over_step);
y_ref = ref(:, 1:T_sim)';
rel_error= zeros(size(y_ref));
% Relative tracking error [%]
for i = 1:length(y_ref)
    for j = 1:ny
        if y_ref(i,j) ~= 0
            rel_error(i,j) = abs(y_sim2(i,j) - y_ref(i,j)) / y_ref(i,j) * 100;
        else
            % Guard against division by zero when the reference is zero
            rel_error(i,j) = 0;
        end
    end
end
Avg_rel_error = mean(rel_error);
for j=1:ny
    RMSE(j) = 1/1 * rmse(ref(j,1:T_sim), y_sim_meas2(1:T_sim,j)');
    % RMSE2(j) = 1/1 * rmse(ref(j,201:T_sim), y_sim_meas2(201:T_sim,j)');
end
disp(RMSE)
% disp(RMSE2)
%% Graphic evaluation
figure
plot(0:T_sim-1, y_sim2, linewidth=1.5)
hold on
plot(0:T_sim-1, ref(:,1:T_sim)', 'g--', linewidth=1)
xlabel('Timestep [-]')
ylabel('Height [cm]')
title("Output", 'FontSize', 14)
legend("DeePC", "reference");
set(gcf, 'Color', 'w');
grid on
exportgraphics(gcf, 'DeePC_general.png', ...
    'Resolution', 600);


figure
subplot(3, 1, 1)
plot(0:T_sim-1, y_sim_meas2, linewidth=1.5)
hold on
plot(0:T_sim-1, ref(:,1:T_sim)', 'g--', linewidth=1)
xlabel('Discrete timestep')
title("Output", 'FontSize', 14)
legend("y", "y_{ref}");
grid on

subplot(3, 1, 2)
plot(0:T_sim-1, u_sim2, 'g')
xlabel('Discrete timestep')
title("Control input", 'FontSize', 14)
grid on

subplot(3, 1, 3)
plot(0:T_sim-1, rel_error, 'r')
xlabel('Discrete timestep')
title("Relative Error", 'FontSize', 14)
grid on



figure
subplot(2,1,1)
plot(0:T_sim-1, y_sim_meas2(:,1), linewidth=1.5)
hold on
plot(0:T_sim-1, y_sim_meas(:,1), linewidth=1.5)
plot(0:T_sim-1, ref(1,1:T_sim)', 'g--', linewidth=1)
plot(0:T_sim-1, 25*ones(1,T_sim)', 'b--', linewidth=1)
% plot(0:T_sim-1, -25*ones(1,T_sim)', 'g--', linewidth=1)
% plot(0:T_sim-1, 15*ones(1,T_sim)', 'g--', linewidth=1)
grid on
title("Velocity (y_1)", 'FontSize', 14)
subplot(2,1,2)
plot(0:T_sim-1, y_sim_meas2(:,2), linewidth=1.5)
hold on
plot(0:T_sim-1, y_sim_meas(:,2), linewidth=1.5)
plot(0:T_sim-1, ref(2,1:T_sim)', 'g--', linewidth=1)
plot(0:T_sim-1, -15*ones(1,T_sim)', 'r--', linewidth=1)
xlabel('Discrete timestep')
title("Climb rate (y_2)", 'FontSize', 14)
legend("Sc-DeePC", "DeePC", "Reference", "Constraint");
grid on




figure('Color','w')

t = tiledlayout(2,1, ...
    'TileSpacing','compact', ...
    'Padding','compact');

% ---------- Subplot 1 ----------
ax1 = nexttile;
hold on
plot(0:T_sim-1, y_sim_meas2(:,1), 'LineWidth', 1.5)
plot(0:T_sim-1, y_sim_meas(:,1), 'LineWidth', 1.5)
plot(0:T_sim-1, ref(1,1:T_sim)', 'g--', 'LineWidth', 1)
plot(0:T_sim-1, 25*ones(1,T_sim)', 'k--', 'LineWidth', 1)

grid on
set(gca, 'FontSize', 12, 'Color', 'w')

ylabel('$y_1\ [\mathrm{ft/s}]$', ...
    'Interpreter', 'latex', 'FontSize', 14)

% ---------- Subplot 2 ----------
ax2 = nexttile;
hold on
plot(0:T_sim-1, y_sim_meas2(:,2), 'LineWidth', 1.5)
plot(0:T_sim-1, y_sim_meas(:,2), 'LineWidth', 1.5)
plot(0:T_sim-1, ref(2,1:T_sim)', 'g--', 'LineWidth', 1)
plot(0:T_sim-1, -15*ones(1,T_sim)', 'k--', 'LineWidth', 1)

grid on
set(gca, 'FontSize', 12, 'Color', 'w')

xlabel('$\mathrm{k} [-]$', ...
    'Interpreter', 'latex', 'FontSize', 14)
ylabel('$y_2\ [\mathrm{ft/s}]$', ...
    'Interpreter', 'latex', 'FontSize', 14)

% ---------- Legend (global, horizontal) ----------
lgd = legend(ax2, ...
    {'$\mathrm{Sc\mbox{-}DeePC}$', ...
     '$\mathrm{DeePC}$', ...
     '$\mathrm{Reference}$', ...
     '$\mathrm{Output\ constraint}$'}, ...
    'Interpreter', 'latex', ...
    'FontSize', 12, ...
    'Orientation', 'horizontal');

% Vector export (PDF) for the paper
exportgraphics(gcf, 'figure.pdf', 'ContentType', 'vector');

% High-resolution raster export (PNG, 600 DPI)
exportgraphics(gcf, 'figure.png', 'Resolution', 600);