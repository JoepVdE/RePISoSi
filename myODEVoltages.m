% -------------------------------------------------------------------------
%  myODEVoltages.m  --  Recover the per-node voltages from the same
%                        linearised system as myODE, after the ODE step.
%  Part of RePISoSi - https://github.com/JoepVdE/RePISoSi  -  License: MIT
%  Author : J.L. Van den Eijnden, 2026
%           Original implementation by M. Mentink (Oct. 2022),
%           extended by J.L. Van den Eijnden and A. Vaskuri (Oct. 2023).
%
%  Inputs : t (s), IArray (A), s (struct, see myODE.m).
%  Output : VGroundVector (V) - voltage at each non-reference node.
% -------------------------------------------------------------------------
function VGroundVector = myODEVoltages(t, IArray, s)   %#ok<INUSL>
    s.signArray = sign(IArray(1:s.numLinesAlongWire));                                          % Direction of the current 
	%s.INormalArray = abs(IArray(1:s.numLinesAlongWire))*ones(1,s.numThermalSubdivisions) - s.ISCArray;  
    s.INormalArray = abs(IArray(1:s.numLinesAlongWire))*ones(1, s.numThermalSubdivisions) - s.IcArray; %Critical state model

	s.INormalArray = s.INormalArray.*(s.INormalArray > 0);                                      % Current travelling in a stabilizer (A)
	
    voltageVectorLongitudinal = -s.signArray.*sum(s.INormalArray.*s.RNormalArray, 2);           % Voltage due to longitudial resistance of a stabilizer (ohm) 
    voltageVectorTrans = -IArray(s.numLinesAlongWire+1:s.numLines).*s.RArrayTrans;              % Voltage due to turn-to-turn resistance (ohm) 
    s.bVector = [voltageVectorLongitudinal; voltageVectorTrans; s.externalVoltagedIdtVector];   
    s.vVector = s.MLPInv*s.bVector;                                                             
    VGroundVector = s.vVector(s.numLines+1:s.numLines+s.numPoints-1);                           % Voltage (V)
end
