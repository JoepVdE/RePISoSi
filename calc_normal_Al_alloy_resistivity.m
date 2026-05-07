% -------------------------------------------------------------------------
%  calc_normal_Al_alloy_resistivity.m  --  Coarse linear approximation of
%                                          aluminium alloy resistivity.
%  Part of RePISoSi - https://github.com/JoepVdE/RePISoSi  -  License: MIT
%  Author : J.L. Van den Eijnden, 2026
%           Original implementation by M. Mentink (Oct. 2022),
%           extended by J.L. Van den Eijnden and A. Vaskuri (Oct. 2023).
%
%  rho_normal = rho0 + rhoT * T  with rho0 = 1e-8 Ohm*m, rhoT = 1e-9.
%  Used as a fallback when the BlochGrun-based table is not appropriate.
%
%  Input : temperatureArray (K), any shape.
%  Output: rho_normal (Ohm*m), same shape.
% -------------------------------------------------------------------------
function rho_normal = calc_normal_Al_alloy_resistivity(temperatureArray)
    rho0 = 1e-8; % Saturates at 1e-8 ohm*m
    rhoT = 1e-9; % Rises at higher temperatures 
    rho_normal = temperatureArray.*rhoT + rho0; 
end