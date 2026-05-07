% -------------------------------------------------------------------------
%  calc_normal_Al10SiMg_resistivity.m  --  AlSi10Mg resistivity from a
%                                          Bloch-Grueneisen fit.
%  Part of RePISoSi - https://github.com/JoepVdE/RePISoSi  -  License: MIT
%  Author : J.L. Van den Eijnden, 2026
%           Original implementation by M. Mentink (Oct. 2022),
%           extended by J.L. Van den Eijnden and A. Vaskuri (Oct. 2023).
%
%  Input : temperatureArray (K), any shape.
%  Output: rho_normal (Ohm*m), same shape.
% -------------------------------------------------------------------------
function rho_normal = calc_normal_Al10SiMg_resistivity(temperatureArray)
    rho_normal = 1e-6 * BlochGrun(temperatureArray, 0.2594, 484, 0.026, 5);
end