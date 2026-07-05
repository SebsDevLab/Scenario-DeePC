%% File Name: system_boundaries_multi.m
% Author: Sebastian Zieglmeier
% Date last updated: 10.06.2026
% Description: Multi-output version of system_boundaries. Checks each element
% of a system variable against its own (vector-valued) lower/upper bound,
% clips any violations, and returns the number of violations in w (rather than
% a 0/1 flag). Used for multi-output systems such as the Boeing 747 benchmark.
% Sources:
%
%
% Inputs:
%   var: Given system variable
%   constraints: The system boundaries
%   var_str: A string to fetch the correct system boundaries
%
% Outputs:
%   var: Given system variable (clipped or original value)
%   w: number of elements that were clipped beyond the threshold
%   warn: warning message for clipping
%
% Notes:
% An overwriting of the state is very common due to physical limitations,
% while the overwriting of the input or the output should not happen that 
% often and can therefore indicate a mistake as the predictive control uses
% the input and output constraints itself. Therefore should u and y be
% normally within the system boundaries


function [var, w, warn] = system_boundaries_multi(var, constraints, var_str)
    var_min = eval(["constraints." + var_str + "_min"]);
    var_max = eval(["constraints." + var_str + "_max"]);
    warn = ["Warning: overwrite due to " + var_str + " outside of system boundary"];
    w = 0; % violation counter
    for i = 1:1:size(var,2)
        if var(1, i) < var_min(i) || var(1, i) > var_max(i)
            var_old = var(1,i);
            var(1, i) = min(max(var_min(i), var(1, i)), var_max(i));
            if abs(var_old - var(1,i)) > 1e-3
                w = w+1; 
                % Only produce warning bool, if the difference is over a 
                % certain threshold, as minimal numerical inaccuracies 
                % in the optimization problem can occur
            end
        end
    end
end