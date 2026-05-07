% -------------------------------------------------------------------------
%  n_value.m  --  Power-law n-value of HTS tape vs. temperature.
%  Part of RePISoSi - https://github.com/JoepVdE/RePISoSi  -  License: MIT
%  Author : J.L. Van den Eijnden, 2026
%           Original implementation by M. Mentink (Oct. 2022),
%           extended by J.L. Van den Eijnden and A. Vaskuri (Oct. 2023).
%
%  Inputs:
%      T           : temperature (K), scalar / vector / matrix
%      HTStapeName : string, currently only Fujikura tapes are supported.
%  Output:
%      n           : n-value (dimensionless), same shape as T.
%                    Saturates at 66 below 4.2 K and at 7 above 77 K.
%
%  TODO: extend to non-Fujikura tapes. The original code used
%        "HTStapeName == 'Fujikura'" which performs an element-wise
%        comparison and silently fails when HTStapeName is longer than 8
%        characters. Replaced with contains() below.
% -------------------------------------------------------------------------
function [n] = n_value(T, HTStapeName)
    if contains(HTStapeName, 'Fujikura')
        n_meas = [7,  11, 12, 13, 15, 20, 34, 66]';
        T_meas = [77, 65, 50, 40, 30, 20, 10, 4.2]';
        n      = interp1(T_meas, n_meas, T);
        n(T < min(T_meas)) = 66;
        n(T > max(T_meas)) = 7;
    else
        % TODO: parameterize for other tape types.
        error('n_value:unsupportedTape', ...
              'n-value table not implemented for tape "%s".', HTStapeName);
    end
end