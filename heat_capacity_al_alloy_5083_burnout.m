% Reference: 
%   LNG Materials and Fluids. Ed. Douglas Mann National Bureau of Standards,
%   Cryogenics Division First Edition, 1977. 

function [C] = heat_capacity_al_alloy_5083_burnout(T)

    rho = 2650; % Density of Al 5083 (kg/m^3)
    m = 2; % Mass (kg) 
    
    %cp = zeros(size(T)); % Init specific heat (J/(kg*K)) 

    %idx = find(T <= 4.0); 
    %cp(idx) = NaN; 

    idx = find(T >= 4.0 & T <= 16.0); 
    cp(idx) = 0.153508701 - 0.104679998*T(idx).^1 + 0.0454540137*T(idx).^2 - 0.00328604666*T(idx).^3 + 1.21017723E-4*T(idx).^4; 

    idx = find(T > 16.0 & T <= 47.0); 
    cp(idx) = -7.85853106 + 1.97771153*T(idx).^1 - 0.1627576*T(idx).^2 + 0.00632458042*T(idx).^3 - 5.22112155E-5*T(idx).^4;

    idx = find(T > 47.0 & T <= 130.0); 
    cp(idx) = -74.8892486 - 1.13922203*T(idx).^1 + 0.176562541*T(idx).^2 - 0.00147785574*T(idx).^3 + 3.9391632E-6*T(idx).^4;

    idx = find(T > 130.0 & T <= 300.0); 
    cp(idx) = -263.527985 + 9.29321629*T(idx).^1 - 0.0121044208*T(idx).^2 - 6.74100976E-5*T(idx).^3 + 1.6547628E-7*T(idx).^4;

    idx = find(T > 300.0); 
    cp(idx) = NaN; 
    
    C = cp.*m; % Heat capacity of support cylinder (J/K) 
    
end
