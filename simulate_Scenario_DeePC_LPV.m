%% File Name: simulate_Scenario_DeePC_LPV.m
% Author: Sebastian Zieglmeier
% Date last updated: 10.06.2026
% Description: Closed-loop simulation of the nonlinear two-tank system
% (LPV model) controlled by DeePC and Scenario-DeePC for a given reference
% trajectory. The script first runs standard DeePC to collect closed-loop
% prediction errors, builds the scenario buffer from those errors, and then
% runs Scenario-DeePC for comparison.
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
sys_name = "LPV_2_Tank";
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
T_sim = 1200;       % Simulation length [steps]
T_fut = 5;          % Prediction horizon N
T_d = 200;          % Number of data points used to build the Hankel matrices
T_ini = 20;         % Length of the initial (past) window T_ini

% DeePC cost-function weights and regularization parameters
lambda_ini = 1e7;   % Penalty on the initial-condition slack sigma_y
lambda_g = 1e4;     % Penalty on the regularizer ||g||^2
r = 1e-2;           % Input weight R
q = 1e4;            % Output (tracking) weight Q

%% Load system
% LPV state-space model used to generate the data for the data-driven control
A_theta = sys.model.A;
B_c = sys.model.B;
C_c = sys.model.C;
D_c = sys.model.D;

A_sim_theta_c = sys.sim_model_1.A_sim;
B_sim_c = sys.sim_model_1.B_sim;
C_sim_c = sys.sim_model_1.C_sim;
D_sim_c = sys.sim_model_1.D_sim;


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
% sequences injected during data collection. They are NOT the paper's
% prediction-error scenarios w^{(j)} (those are built later from y_err).
sigma_meas = 0.01 * sys.constraints.y_max;  % measurement-noise std, scaled to max output
sigma_proc = 0.1*sigma_meas;                % process-noise std
rng(2,'twister');                           % fix the RNG seed for reproducibility
scaling_factor = 1;
for i = 1:T_d
    u_data_sys(i) = rand(1)*scaling_factor*i;               % persistently exciting excitation input
    if u_data_sys(i) >= constraints.u_max
        u_data_sys(i) = mod(u_data_sys(i), constraints.u_max);
        if u_data_sys(i) < constraints.u_max/2
            u_data_sys(i) = u_data_sys(i) + constraints.u_max/2;
        end
    end
    [u_data_sys(i, :), w, warn] = system_boundaries(u_data_sys(i, :), sys.constraints, "u");
    if w == 1
        disp(warn);
    end

    if active_noise == 2
        v_meas(i,:) = sigma_meas * randn(1, ny);
        v_proc(i,:) = sigma_proc * randn(1, nx);
    elseif active_noise == 3
        v_meas(i,:) = sigma_meas * randn(1, ny);
        v_proc(i,:) = zeros(1, nx);
    elseif active_noise == 4
        v_meas(i,:) = zeros(1, ny);
        v_proc(i,:) = sigma_proc * randn(1, nx);
    else
        v_meas(i,:) = zeros(1, ny);
        v_proc(i,:) = zeros(1, nx);
    end

    [A, B, C, D] = discretize_LPV(A_theta, B_c, C_c, D_c, x_data_sys(i, :), sys.T_samp);

    x_data_sys(i+1, :) = (A * x_data_sys(i, :)' + B * u_data_sys(i))' + v_proc(i,:); 
    y_data_sys(i,:) = C * x_data_sys(i, :)' + D * u_data_sys(i); 
    y_meas(i,:) = y_data_sys(i,:) + v_meas(i,:);
    
    % Enforce state and output constraints on the collected data
    [x_data_sys(i+1, :), w, warn] = system_boundaries(x_data_sys(i+1, :), sys.constraints, "x");
    if w == 1
        disp(warn);
    end
    [y_data_sys(i, :), w, warn] = system_boundaries(y_data_sys(i, :), sys.constraints, "y");
    if w == 1
        disp(warn);
    end
    % Collecting data in a certain (approx.) range via scaling_factor
    if y_data_sys(i,:) < 10
        scaling_factor = 1;
    elseif y_data_sys(i,:) > 20 
        scaling_factor = .1;
    end
end
u_data = u_data_sys;
y_data = y_meas;

%% Diagnostics: collected data and measurement-noise distribution
figure
plot(y_data_sys);
hold on;
plot(y_data);


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
num_hankel_cols = T_d - L;
H_u = u_data(1:L);
H_y = y_data(1:L);
for i = 2:num_hankel_cols+1
    H_u = [H_u, u_data(i:i+L-1)];
    H_y = [H_y, y_data(i:i+L-1)];
end





%% Initialize controller and simulation variables (standard DeePC)

if control_name == "MPC"
    control = MPC(T_d, T_ini, T_fut, r, q, lambda_ini, lambda_g, H_u, H_y, sys);
elseif control_name == "DeePC_fast"
    control = DeePC_fast(T_d, T_ini, T_fut, r, q, lambda_ini, lambda_g, H_u, H_y, sys);
else
    disp("Warning: No controller selected. Check spelling of control_name.");
end

% Preallocate arrays to store simulation results
u_sim = zeros(T_sim, 1);
x_sim = zeros(T_sim + 1, nx);
y_sim = zeros(T_sim, 1);
y_sim_pred = zeros(T_sim, ny);
y_sim_meas = zeros(T_sim, ny);
w_sim_meas = zeros(T_sim, ny);
w_sim_proc = zeros(T_sim, nx);

y_D = zeros(T_sim, 1);

%% Initial condition:
x_sim(1, :) = zeros(nx, 1);
u_past_sim = zeros(T_ini, 1);
y_past_sim = zeros(T_ini, 1);

%% Reference Trajectory

ref_name = "Smooth_Step";
ref = get_ref2(ref_name, T_sim, ny, T_fut, T_ini,15,15);
% Multi-level staircase reference (three operating points: 15 -> 20 -> 25 cm)
ref(290:310) = linspace(15,20,21);
ref(311:590) = 20;
ref(590:610) = linspace(20,25,21);
ref(611:end)= 25;
% Available reference profiles: Smooth_Step_to_Sinus, Smooth_Step, Sinus, Step
% Reference values can be changed in get_ref.m
over_step = 0;
%% Simulation loop (standard DeePC)
for i = 1:T_sim 
    y_reference = ref(i:i+T_fut-1);
    [u_fut, y_fut, g_fut, sigma_fut, opt_error_flag] = control.step(u_past_sim, y_past_sim, y_reference', H_u, H_y); 
    if opt_error_flag == 1
        disp(["Warning: Opt_error, Timestep: " + string(i)]);
    end
    u_sim(i) = u_fut(1); % Apply only the first input (receding horizon)
    % Enforce input constraints
    [u_sim(i, :), w, warn] = system_boundaries(u_sim(i, :), sys.constraints, "u");
    if w == 1
        disp(warn);
    end

    if active_noise == 2
        w_sim_meas(i,:) = sigma_meas * randn(1, ny);
        w_sim_proc(i,:) = sigma_proc * randn(1, nx);
    elseif active_noise == 3
        w_sim_meas(i,:) = sigma_meas * randn(1, ny);
        w_sim_proc(i,:) = zeros(1, nx);
    elseif active_noise == 4
        w_sim_meas(i,:) = zeros(1, ny);
        w_sim_proc(i,:) = sigma_proc * randn(1, nx);
    else
        w_sim_meas(i,:) = zeros(1, ny);
        w_sim_proc(i,:) = zeros(1, nx);
    end


    [A_sim, B_sim, C_sim, D_sim] = discretize_LPV(A_sim_theta_c, B_sim_c, C_sim_c, D_sim_c, x_sim(i, :), sys.T_samp);

    % Simulate the real (LPV) system
    x_sim(i+1, :) = (A_sim * x_sim(i, :)' + B_sim * u_sim(i))' + w_sim_proc(i,:);
    y_sim(i) = C_sim * x_sim(i, :)' + D_sim * u_sim(i);
    y_sim_pred(i) = y_fut(1,:);
    y_sim_meas(i) = y_sim(i) + w_sim_meas(i,:);
    
    % Enforce state constraints
    [x_sim(i+1, :), w, warn] = system_boundaries(x_sim(i+1, :), sys.constraints, "x"); 
    if w == 1
        disp(warn);
    end
    % Count output-constraint violations (do not clip the measured output)
    [~, w, warn] = system_boundaries(y_sim_meas(i, :), sys.constraints, "y"); 
    if w == 1
        over_step = over_step + 1;
    end

    % Store the measured output for the past window
    y_D(i) = y_sim_meas(i);
    % Update the past input/output windows for the next step
    u_past_sim = [u_past_sim(2:end); u_sim(i)];
    y_past_sim = [y_past_sim(2:end); y_D(i)];
end
disp("Constraint violations without Scenario approach: ");
disp(over_step);
%% Diagnostics: empirical distributions of noise and prediction errors (DeePC)
% Histograms with a fitted Gaussian for: injected measurement noise (w_sim_meas),
% the noise-free prediction error e_1, and the measured prediction error e_2.
% e_2 is reused below to build the scenario buffer, so this block must run.
figure
subplot(3,1,1)
histogram(w_sim_meas, 'Normalization', 'pdf');
hold on;
w=w_sim_meas;
mu_hat = mean(w);
sigma_hat = std(w);
xline(mu_hat, "r--", 'LineWidth', 2);
% x-grid for the fitted Gaussian
x_vals = linspace(min(w), max(w), 200);

% Gaussian PDF
pdf_vals = (1/(sigma_hat*sqrt(2*pi))) * ...
           exp(-0.5*((x_vals - mu_hat)/sigma_hat).^2);
% Plot fitted Gaussian
plot(x_vals, pdf_vals, 'LineWidth', 2);

e_1 = y_sim - y_sim_pred;             % noise-free prediction error
subplot(3,1,2)
histogram(e_1, 'Normalization', 'pdf');
hold on;
w=e_1;
mu_hat = mean(w);
sigma_hat = std(w);
xline(mu_hat, "r--", 'LineWidth', 2);
% x-grid for the fitted Gaussian
x_vals = linspace(min(w), max(w), 200);

% Gaussian PDF
pdf_vals = (1/(sigma_hat*sqrt(2*pi))) * ...
           exp(-0.5*((x_vals - mu_hat)/sigma_hat).^2);

% Plot fitted Gaussian
plot(x_vals, pdf_vals, 'LineWidth', 2);
e_2 = y_sim_meas - y_sim_pred;        % measured prediction error (used for scenarios)
subplot(3,1,3)
histogram(e_2, 'Normalization', 'pdf');
hold on;
w=e_2;
mu_hat = mean(w);
sigma_hat = std(w);
xline(mu_hat, "r--", 'LineWidth', 2);
% x-grid for the fitted Gaussian
x_vals = linspace(min(w), max(w), 200);

% Gaussian PDF
pdf_vals = (1/(sigma_hat*sqrt(2*pi))) * ...
           exp(-0.5*((x_vals - mu_hat)/sigma_hat).^2);

% Plot fitted Gaussian
plot(x_vals, pdf_vals, 'LineWidth', 2);

%% Numerical evaluation (DeePC)
y_ref = ref(:, 1:T_sim)';
rel_error= zeros(size(y_ref));
% Relative tracking error [%]
for i = 1:length(y_ref)
    if y_ref(i) ~= 0
        rel_error(i) = abs(y_sim_meas(i) - y_ref(i)) / y_ref(i) * 100;
    else
        % Guard against division by zero when the reference is zero
        rel_error(i) = 0;
    end
end
Avg_rel_error = mean(rel_error);
for j=1:ny
    RMSE(j) = 1/1 * rmse(ref(j,1:200), y_sim_meas(1:200,j)');
    RMSE2(j) = 1/1 * rmse(ref(j,201:1000), y_sim_meas(201:1000,j)');
end
disp(RMSE)
disp(RMSE2)
%% Graphic evaluation (DeePC)
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

%% Per-operating-point ACF of the two-tank prediction error (buffer stride M)
M      = 2;                    % Correlation horizon: buffer keeps every (M+1)-th error
segs   = [  50  250;
           350  550;
           650 1150];
w_err  = e_2;                  % Measured prediction error w^{(j)} (T_sim x 1)
nSeg   = size(segs,1);

% Thin each OP segment by (M+1), then de-mean to remove that OP's bias
seg_w = cell(nSeg,1);
for sgi = 1:nSeg
    s = w_err(segs(sgi,1):segs(sgi,2));
    s = s(1:M+1:end);                       % stride M+1
    seg_w{sgi} = s - mean(s);
end
N_op   = cellfun(@numel, seg_w);
maxLag = min(30, floor(min(N_op)/4));        % keep lags well-supported after thinning

% per-OP ACF of the thinned stream
acf_op = zeros(maxLag+1, nSeg);
for sgi = 1:nSeg
    wc = seg_w{sgi};
    d  = sum(wc.^2);
    for tau = 0:maxLag
        acf_op(tau+1,sgi) = sum(wc(1:end-tau).*wc(1+tau:end)) / d;
    end
end

% ---------- paper figure: one subplot per operating point ----------
lags  = (0:maxLag)';
cBlue = [0 0.447 0.741];
cGray = [0.3 0.3 0.3];

fig = figure('Color','w','Units','centimeters','Position',[2 2 17 13]);
tl  = tiledlayout(nSeg,1,'TileSpacing','compact','Padding','compact');

for sgi = 1:nSeg
    ax = nexttile; hold(ax,'on');
    a    = acf_op(:,sgi);
    conf = 1.96 / sqrt(N_op(sgi));

    plot([0 maxLag],[0 0],'-','Color',[0.5 0.5 0.5],'LineWidth',0.5);
    hB = plot([0 maxLag],[ conf  conf],'--','Color',cGray,'LineWidth',1.0);
         plot([0 maxLag],[-conf -conf],'--','Color',cGray,'LineWidth',1.0);
    stem(lags, a,'filled','Color',cBlue,'MarkerFaceColor',cBlue, ...
         'MarkerSize',4,'LineWidth',1.0,'BaseValue',0);

    set(ax,'FontSize',12,'LineWidth',0.9,'Box','on','TickLabelInterpreter','latex', ...
           'XGrid','on','YGrid','on','GridAlpha',0.12);
    xlim([-0.5 maxLag+0.5]); ylim([min(-0.3,min(a)-0.05) 1.05]); xticks(0:5:maxLag);
    ylabel('$\rho(\tau)$','Interpreter','latex','FontSize',14);
    title(sprintf('OP %d: $k\\in[%d,\\,%d]$ ($N=%d$)',sgi,segs(sgi,1),segs(sgi,2),N_op(sgi)), ...
          'Interpreter','latex','FontSize',12);

    if sgi == 1
        legend(hB,'95\% band','Interpreter','latex', ...
               'Location','northeast','Box','off','FontSize',10);
    end
    if sgi == nSeg
        xlabel('Lag $\tau$','Interpreter','latex','FontSize',14);
    end
end
title(tl,sprintf('Prediction-error autocorrelation per operating point ($M=%d$)',M), ...
      'Interpreter','latex','FontSize',14);


% -------------------------------------------------------------------------
%% Scenario-DeePC run
% -------------------------------------------------------------------------
clear w_scen w_flat
%% Initialize controller and simulation variables (Scenario-DeePC)
N_scen = 20;        % Number of scenarios N_scen
T_fut = 5;          % Prediction horizon N
T_d = 200;          % Number of data points used to build the Hankel matrices
T_ini = 20;         % Length of the initial (past) window T_ini

% DeePC cost-function weights and regularization parameters
lambda_ini = 1e7;   % Penalty on the initial-condition slack sigma_y
lambda_g = 1e4;     % Penalty on the regularizer ||g||^2
r = 1e-2;           % Input weight R
q = 1e4;            % Output (tracking) weight Q

% Last argument 1e6 is the exact-penalty weight on the scenario slack
control = Scenario_DeePC_multi_Sc_cost(T_d, T_ini, T_fut, r, q, lambda_ini, lambda_g, 1e6, H_u, H_y, sys, N_scen);


% Preallocate arrays to store simulation results
u_sim2 = zeros(T_sim, 1);
x_sim2 = zeros(T_sim + 1, nx);
y_sim2 = zeros(T_sim, 1);
y_sim_pred2 = zeros(T_sim, 1);
y_sim_meas2 = zeros(T_sim, 1);
w_sim_meas2 = w_sim_meas;       % reuse identical noise realization for a fair comparison
w_sim_proc2 = w_sim_proc;       % reuse identical noise realization for a fair comparison
error_buffer = zeros(2*N_scen, ny);   % scenario buffer E_buffer of one-step errors w^{(j)}
y_D2 = zeros(T_sim, 1);

%% Initial condition
x_sim2(1, :) = zeros(nx, 1);
u_past_sim = zeros(T_ini, 1);
y_past_sim = zeros(T_ini, 1);

%% Reference trajectory (reused from the DeePC run)
ref_name = "Smooth_Step";
% Available reference profiles: Smooth_Step_to_Sinus, Smooth_Step, Sinus, Step
% Reference values can be changed in get_ref.m
disp("Scenario approach start");
over_step = 0;
%% Simulation loop
for i = 1:T_sim 
    y_reference = ref(i:i+T_fut-1);
    % Sample N_scen scenarios over the horizon from the buffer E_buffer
    idx = randi(size(error_buffer, 1), T_fut * N_scen, 1);
    w_flat = error_buffer(idx, :);
    w_scen = reshape(w_flat, T_fut, ny, N_scen);
    [u_fut, y_fut, g_fut, sigma_fut, opt_error_flag] = control.step(u_past_sim, y_past_sim, y_reference', H_u, H_y, w_scen); 
    if opt_error_flag == 1
        disp(["Warning: Opt_error, Timestep: " + string(i)]);
    end
    u_sim2(i) = u_fut(1); % Apply only the first input (receding horizon)
    % Enforce input constraints
    [u_sim2(i, :), w, warn] = system_boundaries(u_sim2(i, :), sys.constraints, "u");
    if w == 1
        disp(warn);
    end

    % Noise is intentionally NOT regenerated here: w_sim_meas2/w_sim_proc2 were
    % set to the same realizations used in the DeePC run above, so both
    % controllers see identical noise and the comparison is fair.

    % Simulate the real (LPV) system
    [A_sim, B_sim, C_sim, D_sim] = discretize_LPV(A_sim_theta_c, B_sim_c, C_sim_c, D_sim_c, x_sim2(i, :), sys.T_samp);

    x_sim2(i+1, :) = (A_sim * x_sim2(i, :)' + B_sim * u_sim2(i))' + w_sim_proc2(i,:);
    y_sim2(i) = C_sim * x_sim2(i, :)' + D_sim * u_sim2(i);
    y_sim_pred2(i) = y_fut(1,:);
    y_sim_meas2(i) = y_sim2(i) + w_sim_meas2(i,:);
    
    % Enforce state constraints
    [x_sim2(i+1, :), w, warn] = system_boundaries(x_sim2(i+1, :), sys.constraints, "x"); 
    if w == 1
        disp(warn);
    end
    % Count output-constraint violations (do not clip the measured output)
    [~, w, warn] = system_boundaries(y_sim_meas2(i, :), sys.constraints, "y"); 
    if w == 1
        over_step = over_step + 1;
    end

    % Store the measured output for the past window
    y_D2(i) = y_sim_meas2(i);
    % Update the past input/output windows for the next step
    u_past_sim = [u_past_sim(2:end); u_sim2(i)];
    y_past_sim = [y_past_sim(2:end); y_D2(i)];
    y_err(i,:) = y_sim_meas2(i,:) - y_sim_pred2(i,:);   % one-step prediction error w^{(j)}
    if i > 50 && mod(i - 51, M+1) == 0      % keep every (M+1)-th error so buffer entries are decorrelated
        error_buffer = [error_buffer(2:end, :); y_err(i,:)];
    end
end
disp("Constraint violations with Scenario approach: ");
disp(over_step);

%% Numerical evaluation (Scenario-DeePC)
y_ref = ref(:, 1:T_sim)';
rel_error2= zeros(size(y_ref));
% Relative tracking error [%]
for i = 1:length(y_ref)
    if y_ref(i) ~= 0
        rel_error2(i) = abs(y_sim_meas2(i) - y_ref(i)) / y_ref(i) * 100;
    else
        % Guard against division by zero when the reference is zero
        rel_error2(i) = 0;
    end
end
Avg_rel_error = mean(rel_error);

for j=1:ny
    RMSE(j) = 1/1 * rmse(ref(j,1:200), y_sim_meas2(1:200,j)');
    RMSE2(j) = 1/1 * rmse(ref(j,201:T_sim), y_sim_meas2(201:T_sim,j)');
end
disp(RMSE)
disp(RMSE2)

%% Graphic evaluation (Scenario-DeePC)
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
plot(0:T_sim-1, y_sim_meas2, linewidth=1.5)
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


%% Comparison figure: DeePC vs Scenario-DeePC (output and relative error)
window = 20;        % moving-average window length

rel_err2_ma = movmean(rel_error2(:,1), window);
rel_err_ma  = movmean(rel_error(:,1), window);
y_sim_meas2_ma = movmean(y_sim_meas2(:,1), window);
y_sim_meas_ma  = movmean(y_sim_meas(:,1), window);
T_sim = 1000;
k = 200:T_sim-1;
k2 = 201:T_sim;

% Base colors
c_red  = [1 0 0];
c_blue = [0 0 1];

% Lightened versions (alpha = 0 keeps the color, alpha = 1 is white)
alpha = 0.6;
c_red_light  = (1-alpha)*c_red  + alpha*[1 1 1];
c_blue_light = (1-alpha)*c_blue + alpha*[1 1 1];

figure('Color','w')

t = tiledlayout(2,1, ...
    'TileSpacing','compact', ...
    'Padding','compact');

% ---------- Subplot 1 ----------
ax1 = nexttile;
hold on

plot(k, y_sim_meas2(k2,1), 'Color', c_red_light,  'LineWidth', 0.75)
plot(k, y_sim_meas(k2,1),  'Color', c_blue_light, 'LineWidth', 0.75)
plot(k, y_sim_meas2_ma(k2,1), 'Color', c_red,  'LineWidth', 2)
plot(k, y_sim_meas_ma(k2,1),  'Color', c_blue, 'LineWidth', 2)
plot(k, ref(1,k2)', 'g--', 'LineWidth', 2)
plot(k, 25*ones(1,T_sim-200)', 'k--', 'LineWidth', 1)

grid on
set(gca, 'FontSize', 16, 'Color', 'w')

ylabel('$h\ [\mathrm{cm}]$', ...
    'Interpreter', 'latex', 'FontSize', 16)

% ---------- Subplot 2 ----------
ax2 = nexttile;
hold on


% --- Raw signals (light, in the background) ---
plot(k, rel_error2(k2,1), ...
    'LineWidth', 0.8, 'Color', c_red_light)

plot(k, rel_error(k2,1), ...
    'LineWidth', 0.8, 'Color', c_blue_light)

% --- Moving averages (dominant, in the foreground) ---
plot(k, rel_err2_ma(k2), ...
    'LineWidth', 2, 'Color', c_red)

plot(k, rel_err_ma(k2), ...
    'LineWidth', 2, 'Color', c_blue)

grid on
set(gca, 'FontSize', 16, 'Color', 'w')

xlabel('$\mathrm{k} [-]$', ...
    'Interpreter', 'latex', 'FontSize', 16)

ylabel('$\mathrm{Relative\ error\ [\%]}$', ...
    'Interpreter', 'latex', 'FontSize', 16)

% ---------- Legend (global, horizontal) ----------
lgd2 = legend(ax2, ...
    {'$\mathrm{Sc\mbox{-}DeePC}\!: e_{\mathrm{rel,MA}}$', ...
     '$\mathrm{DeePC}\!: e_{\mathrm{rel,MA}}$', ...
     '$\mathrm{Sc\mbox{-}DeePC}\!: e_{\mathrm{rel}}$', ...
     '$\mathrm{DeePC}\!: e_{\mathrm{rel}}$', ...
     }, ...
    'Interpreter', 'latex', ...
    'FontSize', 14, ...
    'Orientation', 'horizontal', ...
    'NumColumns', 2);
lgd = legend(ax1, ...
    {'$\mathrm{Sc\mbox{-}DeePC}\!: y$', ...
     '$\mathrm{DeePC}\!: y$', ...
     '$\mathrm{Sc\mbox{-}DeePC}\!: y_{\mathrm{MA}}$', ...
     '$\mathrm{DeePC}\!: y_{\mathrm{MA}}$', ...
     '$\mathrm{Reference}$', ...
    '$\mathrm{Output\ constraint}$', ...
    }, ...
    'Interpreter', 'latex', ...
    'FontSize', 14, ...
    'Orientation', 'horizontal', ...
    'NumColumns', 2);

% Vector export (PDF) for the paper
exportgraphics(gcf, 'figure_LPV2.pdf', 'ContentType', 'vector');

% High-resolution raster export (PNG, 600 DPI)
exportgraphics(gcf, '2_Tank_LPV2.png', 'Resolution', 600);


%% Per-operating-point RMSE and constraint-violation counts
disp("Sc-DeePC")
size(find(y_sim_meas2(1:T_sim)>sys.constraints.y_max),1) % number of output-constraint violations
for j=1:ny
    RMSE(j) = 1/1 * rmse(ref(j,201:600), y_sim_meas2(201:600,j)');
    RMSE2(j) = 1/1 * rmse(ref(j,601:T_sim), y_sim_meas2(601:T_sim,j)');
end
disp(RMSE)
disp(RMSE2)
disp("DeePC")
size(find(y_sim_meas(1:T_sim)>sys.constraints.y_max),1) % number of output-constraint violations
for j=1:ny
    RMSE(j) = 1/1 * rmse(ref(j,201:600), y_sim_meas(201:600,j)');
    RMSE2(j) = 1/1 * rmse(ref(j,601:T_sim), y_sim_meas(601:T_sim,j)');
end
disp(RMSE)
disp(RMSE2)