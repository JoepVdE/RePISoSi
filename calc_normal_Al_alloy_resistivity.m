%% ***** Aluminium alloy resistivity vs. temperature ***** 
% 
% Inputs: 
% T = temperature (K), can be a matrix 
% Ic0 = critical current (A) near 0 K, scalar 
% Tc = Critical temperature (K), scalar
% 
% Output: 
% Ic = critical current (A) 
% ***** 


function rho_normal = calc_normal_Al_alloy_resistivity(temperatureArray) 
    rho0 = 1e-8; % Saturates at 1e-8 ohm*m
    rhoT = 1e-9; % Rises at higher temperatures 
    rho_normal = temperatureArray.*rhoT + rho0; 
end