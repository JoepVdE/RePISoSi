% ***** Ordinary Differential Equation solving V ***** 
% 
% Solves linearized system v = inv(A)*b   
% 
% Inputs: 
% t = time (s) 
% IArray = current (A) 
% s = struct with multiple fields 
% 
% Output: 
% dI_per_dt = current rate (A/s) 
% 
% ***** 
function VGroundVector = myODEVoltages(t, IArray, s)
    s.signArray = sign(IArray(1:s.numLinesAlongWire));                                          % Direction of the current 
	s.INormalArray = abs(IArray(1:s.numLinesAlongWire))*ones(1, s.numThermalSubdivisions) - s.IcArray;
	s.INormalArray = s.INormalArray.*(s.INormalArray > 0);                                      % Current travelling in a stabilizer (A)
	
    voltageVectorLongitudinal = -s.signArray.*sum(s.INormalArray.*s.RNormalArray, 2);           % Voltage due to longitudial resistance of a stabilizer (ohm) 
    voltageVectorTrans = -IArray(s.numLinesAlongWire+1:s.numLines).*s.RArrayTrans;              % Voltage due to turn-to-turn resistance (ohm) 
    s.bVector = [voltageVectorLongitudinal; voltageVectorTrans; s.externalVoltagedIdtVector];   
    s.vVector = s.MLPInv*s.bVector;                                                             
    VGroundVector = s.vVector(s.numLines+1:s.numLines+s.numPoints-1);                           % Voltage (V)
end
