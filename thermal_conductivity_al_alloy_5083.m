% -------------------------------------------------------------------------
%  thermal_conductivity_al_alloy_5083.m  --  Thermal conductivity of
%                                            Al-5083 (NBS cryogenic fit).
%  Part of RePISoSi - https://github.com/JoepVdE/RePISoSi  -  License: MIT
%  Author : J.L. Van den Eijnden, 2026
%           Original implementation by M. Mentink (Oct. 2022),
%           extended by J.L. Van den Eijnden and A. Vaskuri (Oct. 2023).
%
%  Reference:
%      LNG Materials and Fluids. Ed. Douglas Mann, National Bureau of
%      Standards, Cryogenics Division, First Edition, 1977.
%
%  k = thermal_conductivity_al_alloy_5083(T) returns W/(m*K) for each
%  element of T (K). Returns NaN for T outside the validated range
%  (currently 4-849 K).
% -------------------------------------------------------------------------
function [k] = thermal_conductivity_al_alloy_5083(T)

    
    k = zeros(size(T)); % Thermal conductivity (W/(m*K)) 

    %idx = find(T <= 4.0); 
    %k(idx) = NaN; 
    
    idx = find(T >= 4 & T <= 130.0); 
    k(idx) = -1.051689 + 0.9981587*T(idx).^1 - 0.004145406*T(idx).^2 + 8.950483E-6*T(idx).^3; 

    idx = find(T > 130.0 & T <= 300.0); 
    k(idx) = 11.817 + 0.6909558*T(idx).^1 - 0.001587888*T(idx).^2 + 1.597797E-6*T(idx).^3; 

    idx = find(T > 300.0 & T <= 849.0); 
    k(idx) = 30.01698 + 0.4662148*T(idx).^1 - 6.570887E-4*T(idx).^2 + 3.182847E-7*T(idx).^3; 

    idx = find(T > 849.0); 
    k(idx) = NaN; 
    
end
