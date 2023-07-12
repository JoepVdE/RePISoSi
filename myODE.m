% ***** Ordinary Differential Equation solving dI/dt = V(t)/L ***** 
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
function dI_per_dt = myODE(t, IArray, s) 
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
