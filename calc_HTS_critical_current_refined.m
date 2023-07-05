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

function Ic = calc_HTS_critical_current_refined(T, Ic0,x,theta_deg,B)
%x = [0.0416305327396655 0.371106284489619 0.0374417653765792 0.0614639207403106 11.7637639979036 30.0684135381674 0.166640442478695 0.235340235784599 1.35754599801058 3.38655866222062 36.0368585670277 0.992352349642806 1.88578369197932 1.63330222942126 -0.00643731012226043 -0.175037084744646 8.11190080396986 1.84577462825066 0.00697856257152212];


    g0 = x(1);      % (dimensionless) 
    g1 = x(2);      % (dimensionless) 
    g2 = x(3);      % (dimensionless) 
    g3 = x(4);      % (dimensionless) 

    Jc0 = x(5);     % (A/mm2)
    Tc = x(6);      % (K) 
    Bc = x(7);      % (T)
    m1 = x(8);      % (dimensionless) 
    n1 = x(9);      % (dimensionless) 

    Jc0ab = x(10);  % (A/mm2)
    Tab = x(11);    % (K) 
    Bab = x(12);    % (T) 
    c = x(13);      % (dimensionless) 
    n2 = x(14);     % (dimensionless) 
    h = x(15);      % (dimensionless) 
    p = x(16);      % (dimensionless) 
    Tab2 = x(17);   % (K) 
    gamma = x(18); 
    offset = x(19); 
    
    theta = theta_deg.*(pi()./180); 
    
    % Anisotropy shape of the critical current density
    g = g0 + g1.*exp(-B.*(g2.*exp(T.*g3))); 

    % Magnetic field dependence of the critical current density in the b direction
    Jcab = Jc0ab.*exp(-(T./Tab).^n2).*((B./Bab) + c).^(h.*(T./Tab2) + p); 

    % Magnetic field dependence of the critical current density
    Jcc = Jc0.*exp(-(T./Tc).^n1).*(exp(-(B./Bc).^m1)); 

    Jc = 1e5.*(offset + real((min(Jcc, Jcab)) + (max(0, Jcab - Jcc)./(1 + ((theta - pi()./2)./g).^(gamma))))); 

    Jc(find(meas_Jc == -1)) = -1;     

    %Ic = Ic0*(Tc - T)./Tc; 
    Ic = Ic.*(Ic > 0);
end