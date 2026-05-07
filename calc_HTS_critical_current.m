% -------------------------------------------------------------------------
%  calc_HTS_critical_current.m  --  Linear T-dependence of HTS Ic.
%  Part of RePISoSi - https://github.com/JoepVdE/RePISoSi  -  License: MIT
%  Author : J.L. Van den Eijnden, 2026
%           Original implementation by M. Mentink (Oct. 2022),
%           extended by J.L. Van den Eijnden and A. Vaskuri (Oct. 2023).
%
%  Ic = calc_HTS_critical_current(T, Ic0, Tc) returns Ic that decreases
%  linearly from Ic0 at 0 K to 0 at Tc. Negative values are clamped to 0.
%
%  Inputs:
%      T   : temperature (K), scalar / vector / matrix
%      Ic0 : critical current near 0 K (A), scalar
%      Tc  : critical temperature (K), scalar
%  Output:
%      Ic  : critical current (A), same shape as T
% -------------------------------------------------------------------------
function Ic = calc_HTS_critical_current(T, Ic0, Tc)
    Ic = Ic0*(Tc - T)./Tc; 
    Ic = Ic.*(Ic > 0);
end