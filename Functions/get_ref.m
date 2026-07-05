%% File Name: get_ref.m
% Author: Sebastian Zieglmeier 
% Date last updated: 10.06.2026
% Description: Builds the reference trajectory for a given profile name.
% Sources:
%
%
% Inputs:
% ref_name: name of the reference profile to generate
% T_sim: Simulation time
% ny: number of outputs 
% T_fut: Prediction horizon for a sufficient number of steps over the simulation horizon
% T_ini: Number values in u_past and y_past of the data-driven component
%
%
% Outputs:
% ref: the reference trajectory
%
% Notes: 
% smooth_len: Smoothening of the step e.g. to satisfy system limitations
% step_len: Length of the Step without overlying sinus 
% sin_period_len: length of the sinus period
% y_0: initial value of reference trajectory
% y_step: final value of the step
% y_sin: amplitude of the sinus

function ref = get_ref(ref_name, T_sim, ny, T_fut, T_ini)
    ini_len = T_ini;
    smooth_len = 15;
    step_len = 20;
    sin_period_len = 40;
    y_0 = 0; 
    y_step = 10; 
    y_sin = 5;
    if ref_name == "Smooth_Step_to_Sinus_paper_LTI"
        ref = Smooth_Step_to_Sinus_paper_LTI(T_sim, ny, T_fut, 19, 6, 21, 40, 0, 10, 5);
    elseif ref_name == "Smooth_Step_to_Sinus_paper_LPV"
        ref = Smooth_Step_to_Sinus_paper_LPV(T_sim, ny, T_fut, 19, 27, 81, 130, 0, 15, 5);
    elseif ref_name == "Smooth_Step_to_Sinus_paper_LPV_robust"
        ref = Smooth_Step_to_Sinus_paper_LPV(T_sim, ny, T_fut, 19, 27, 81, 130, 0, 25, 5);
    else
        ref = eval([ref_name + "(T_sim, ny, T_fut, ini_len, smooth_len, step_len, sin_period_len, y_0, y_step, y_sin)"]);
    end
end