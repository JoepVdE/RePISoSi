% -------------------------------------------------------------------------
%  calc_HTS_critical_current_refined.m  --  *** BROKEN / DEAD CODE ***
%
%  Part of RePISoSi - https://github.com/JoepVdE/RePISoSi  -  License: MIT
%  Author : J.L. Van den Eijnden, 2026
%           Original implementation by M. Mentink (Oct. 2022),
%           extended by J.L. Van den Eijnden and A. Vaskuri (Oct. 2023).
%
%  WARNING: This function references the symbols `meas_Jc` and `Ic0` that
%  are not in its argument list and not defined in scope. It cannot run.
%  It is preserved here only because the body documents the intended
%  multi-parameter Jc(B, T, theta) parametrisation. For a working
%  implementation use parametrisation_fujikura.m instead.
%
%  TODO: rewrite this function (or delete) once the refined parametrisation
%        is properly integrated. Until then it is gated by an early error.
% -------------------------------------------------------------------------
function Ic = calc_HTS_critical_current_refined(T, Ic0, x, theta_deg, B)
    error('calc_HTS_critical_current_refined:notImplemented', ...
          ['calc_HTS_critical_current_refined is broken (references ', ...
           'undefined symbols meas_Jc / Ic0). Use parametrisation_fujikura instead.']);

    %#ok<*UNRCH>  %#ok<*INUSD>
    % --- Original (non-functional) body retained below for reference ---


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