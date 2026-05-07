% -------------------------------------------------------------------------
%  checkconductivity.m  --  One-line smoke test for thermal_conductivity_Al10SiMg.
%  Part of RePISoSi - https://github.com/JoepVdE/RePISoSi  -  License: MIT
%  Author : J.L. Van den Eijnden, 2026
%           Original implementation by M. Mentink (Oct. 2022),
%           extended by J.L. Van den Eijnden and A. Vaskuri (Oct. 2023).
%
%  Evaluates the thermal-conductivity model on a 5x5 patch at T = 5 K and
%  prints the result. Useful for quickly verifying that the BlochGrun fit
%  is loaded correctly.
% -------------------------------------------------------------------------

thermal_conductivity_Al10SiMg(ones(5, 5) * 5)
