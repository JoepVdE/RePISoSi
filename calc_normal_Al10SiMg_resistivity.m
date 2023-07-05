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


function rho_normal = calc_normal_Al10SiMg_resistivity(temperatureArray) 
    rho_normal = 1E-6*BlochGrun(temperatureArray,0.2594,484,0.026,5);
end