% -------------------------------------------------------------------------
%  myODE.m  --  Right-hand side of the line-current ODE system.
%               Solves   dI/dt = inv(L(t)) * V(t)
%               for the line-element currents I.
%  Part of RePISoSi - https://github.com/JoepVdE/RePISoSi  -  License: MIT
%  Author : J.L. Van den Eijnden, 2026
%           Original implementation by M. Mentink (Oct. 2022),
%           extended by J.L. Van den Eijnden and A. Vaskuri (Oct. 2023).
%
%  Inputs:
%      t      : time (s), used by the stiff ODE solver but not in the rhs
%      IArray : current vector (A)  - length numLines
%      s      : struct with fields MLPInv, IcArray, RNormalArray, ...
%  Output:
%      dI_per_dt : dI/dt (A/s) for every line element
% -------------------------------------------------------------------------
function dI_per_dt = myODE(t, IArray, s)   %#ok<INUSL>
    s.signArray = sign(IArray(1:s.numLinesAlongWire));                                          % Direction of the current 
	s.INormalArray = abs(IArray(1:s.numLinesAlongWire))*ones(1, s.numThermalSubdivisions) - s.IcArray; %critical state model
	s.INormalArray = s.INormalArray.*(s.INormalArray > 0);                                      % Current travelling in a normal metal (A)
	
    voltageVectorLongitudinal = -s.signArray.*sum(s.INormalArray.*s.RNormalArray, 2);           % Voltage due to longitudial resistance of a stabilizer (ohm) 
    voltageVectorTrans = -IArray(s.numLinesAlongWire+1:s.numLines).*s.RArrayTrans;              % Voltage due to turn-to-turn resistance (ohm) 
    s.bVector = [voltageVectorLongitudinal; voltageVectorTrans; s.externalVoltagedIdtVector];   % Voltage components  
    s.vVector = s.MLPInv*s.bVector;                                                             % dI/dt = inv(L(t))*V(t)
    s.VGroundVector = s.vVector((s.numLines + 1):(s.numLines + s.numPoints - 1));               % Voltage (V)
    dI_per_dt = s.vVector(1:s.numLines);                                                        % Current rate dI/dt (A/s)
end 
