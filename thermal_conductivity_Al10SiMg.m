% -------------------------------------------------------------------------
%  thermal_conductivity_Al10SiMg.m  --  Thermal conductivity of AlSi10Mg
%                                       via the Wiedemann-Franz law.
%  Part of RePISoSi - https://github.com/JoepVdE/RePISoSi  -  License: MIT
%  Author : J.L. Van den Eijnden, 2026
%           Original implementation by M. Mentink (Oct. 2022),
%           extended by J.L. Van den Eijnden and A. Vaskuri (Oct. 2023).
%
%  k(T) = L * T / rho(T)   with the Lorenz number L = 2.44e-8 V^2/K^2
%  and rho(T) from a Bloch-Grueneisen fit (BlochGrun.m).
% -------------------------------------------------------------------------
function [k] = thermal_conductivity_Al10SiMg(T)
    L = 2.44e-8;   % Lorenz number, V^2 K^-2
    k = 1e6 .* T .* L ./ BlochGrun(T, 0.24233, 432.43717, 0.0247015, 5);
end
