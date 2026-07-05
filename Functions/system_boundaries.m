%% File Name: system_boundaries.m
% Author: Sebastian Zieglmeier 
% Date last updated: 10.06.2026
% Description: Checks whether a system variable satisfies its constraints.
% If a bound is exceeded the value is clipped to the boundary and a warning
% flag is raised.
% Sources: 
%
%
% Inputs:
%   var: Given system variable
%   constraints: The system boundaries
%   var_str: A string to fetch the correct system boundaries
%
% Outputs:
%   var: Given system variable (overwritten or original value)
%   w: boolean variable for overwriting
%   warn: warning for overwriting
%
% Notes: 
% An overwriting of the state is very common due to physical limitations,
% while the overwriting of the input or the output should not happen that 
% often and can therefore indicate a mistake as the predictive control uses
% the input and output constraints itself. Therefore should u and y be
% normally within the system boundaries


function [var, w, warn] = system_boundaries(var, constraints, var_str)
    var_min = eval(["constraints." + var_str + "_min"]);
    var_max = eval(["constraints." + var_str + "_max"]);
    warn = ["Warning: overwrite due to " + var_str + " outside of system boundary"];
    w = 0; % warning bool
    for i = 1:1:size(var,2)
        if var(1, i) < var_min || var(1, i) > var_max
            var_old = var(1,i);
            var(1, i) = min(max(var_min, var(1, i)), var_max);
            if abs(var_old - var(1,i)) > 1e-3
                w = 1; 
                % Only produce warning bool, if the difference is over a 
                % certain threshold, as minimal numerical inaccuracies 
                % in the optimization problem can occur
            end
        end
    end
end