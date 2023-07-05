%% ***** Critical current vs. temperature ***** 
% 
% Inputs: 
% T = temperature (K), can be a matrix 
% Ic0 = critical current (A) near 0 K, scalar 
% Tc = Critical temperature (K), scalar
% 
% Output: 
% Ic = critical current (A) 
% ***** 


function Ic = calc_HTS_critical_current(T, Ic0, Tc)
    Ic = Ic0*(Tc - T)./Tc; 
    Ic = Ic.*(Ic > 0);
end