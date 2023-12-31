% ***** Quench analysis of partially insulated coil *****
% 
%
% Program written by M. Mentink, A. Vaskuri and J.L. Van den Eijnden (Oct., 2023). 
%
% *****


clear all; close all; clc; 

set(0, 'defaultTextFontName', 'times', 'defaultTextFontSize', 14);
set(0, 'defaultAxesFontName', 'times', 'defaultAxesFontSize', 14);

format long
HTStapeName = 'Fujikura FESC-04 4 mm wide'; %name of tape used
I0 = 4431;                    % Driving current (A)  test2AnnaJoep
writedata = 1; %write output file 1 = yes 0 = no
writevideo = 0;

Description = '15turn4431A_80VHeater_allon';
FileName = append(Description,string(datetime('now','TimeZone','local','Format','d-MMM-y')),'_',string(datetime('now','TimeZone','local','Format','HH.mm.ss')),'.txt');  % 

if writedata ==1
fileID = fopen(FileName,'w');
end

if writevideo ==1 
v = VideoWriter(append(Description,'.mp4'),'MPEG-4');
v.FrameRate = 5;
open(v);
end


Helium_cooling = 1; %1 is on, 0 is off.
copper_heat = 1; %1 is on, 0 is off.
copper_conduction = 1;%1 is on, 0 is off.


toffset = [8.071,48.795,70.456,84.040,69.285,100]; % [s] time heater is on
theater = [1.633,5.478,10.304,15.42,19.826]; %[s]   time in between heating pulses
RHeater = 350; % resistance of heater [Ohm]
VHeater = 80; % voltage over heater V

PHeater = VHeater.^2/RHeater;

dTMax = 1;                   % Maximum temperature difference (K) if too large (at 10 K), temperature becomes unaccurate


initial_temperature = 20; % Initial temperature of 10 K
temperature_helium = initial_temperature; %[K]
k_helium = conductivity_helium(temperature_helium);

% ***** Solenoid properties ***** 
RSol = 0.12;                 % Solenoid radius (m)
N_shorts = 19;               % Number of shorts per turn

arclengthSlot = 4*10^(-2);      % arclength of slot (m)
thSlot = 3.7*10^(-3);     %arclength of short (m)

resolution = 4; % number of lines between axial shorts, must be even
resolution = round(resolution/2)*2; %force to be even

numPointsPerTurn = N_shorts*resolution-resolution/2; %make resolution such that shorts are such that they maximalise the axial resistance (off by one turn)

wCond = 0.0085;              % Conductor width (m) 
thReBCO = 2.2E-6;            %thickness ReBCO layer (m)
thCond = 0.0038;             % Effective conductor thickness (including Al stabiliser) (m), 
                             % 3.8 mm = 1.5 mm * 5 mm + 4.5 mm * 3 mm + 1.5 mm * 5 mm)/(1.5 mm + 4.5 mm + 1.5 mm)   

thWallatShort = 5e-3; %[m] thickness of the short wall in between slots. 
heightSlot = 1e-3; %[m] 
ACond = wCond.*thCond;       % Cross-sectional area of the conductor (m)        is actually 7.5 mm but needs to be 8.5 here to include the 1 mm gap

numWindings = 5;            % Number of turn forming a solenoid (dimensionless)
len = numWindings*(wCond + thSlot);     % Length of a solenoid (m)

ignoreMutualInductancesTransverseElements = 0; % (yes = 1, no = 0)
numThermalSubdivisions = 3;  % Number of subdivisions
temperatureSubSteps = 200;   % Number of temperature timesteps per magnetic timestep 




% ***** Current lead properties ******
NRings = 2;
DouterLead = 0.307; %diameter in meter
DinnerLead = 0.208;


ARing = pi*(DouterLead^2-DinnerLead^2)/4; %[m^2]

thRing = 0.01; %[m] thickness copper current lead

crossSectionRing = thRing*(DouterLead-DinnerLead); %[m^2]
rhoCopper = 8960;  %[kg/m^3]

massRing = ARing*thRing*rhoCopper; %[kg]
temperatureRing = initial_temperature; %[K]

RingRRR = 160; %[-] Ekin for annealed OFE copper. Normal copper RRR~110, not ennealed OFE RRR~45


% ***** Superconducting cable properties ***** 
N_tapes = 4;
tapewidth = 4E-3; %[m]
SCthickness = 2.2E-6; %[m]



% ***** Look up tables ***** 
Tinitialise = 0.1:0.1:400;
L = 2.44E-8; %[V^2K^-2]

ktable = 1E6.*Tinitialise.*L./BlochGrun(Tinitialise,0.24233,432.43,0.0247015,5);
rhotable = BlochGrun(Tinitialise,0.24233,432.43717,0.0247015,5)/1E6;
% ***** 
%

% ***** Initialize mesh for FEM ***** 

pointAngle = 2*pi/numPointsPerTurn;  % Angle between the points (rad)

RSolMin = RSol*cos(pointAngle); % Minimum distance from the center due to finite number of points (m)
RSolAvg = 0.5*(RSolMin + RSol); % Average distance (m)



RSolEff = RSol^2/RSolAvg;       % Effective distance (m)


mutualInductanceSpaceRatio = 10; 
mutualInductanceSpaceRatioLimit = 1000;

numPoints = floor(numPointsPerTurn*numWindings + 1);            % Total number of 3D points
numThermalElements = (numPoints - 1)*numThermalSubdivisions;    % Total number of thermal elements

numLinesAlongWire = numPoints - 1; 

dIdt = zeros(numLinesAlongWire,1); %initial dIdt
%numTransverse = floor((numWindings - 1)*numPointsPerTurn + 1); 
numTransverse = floor((numWindings - 1)*N_shorts + 1)-round((numWindings/2)); %It has 0.5 turn leftover each round because shorts spacing is off by one turn. That's why
numTransverse = max(numTransverse, 0); 
numLines = numLinesAlongWire + numTransverse; 

IArrayTrans = zeros(numTransverse-1,1);

% Generate line elements in helix form
[xArray, yArray, zArray] = ... 
    generate_line_elements(numPoints, numPointsPerTurn, RSolEff, numWindings, wCond); 



% Gather space for lines
[x1Array, x2Array, y1Array, y2Array, z1Array, z2Array, pointIndex1, pointIndex2] = ... 
    gather_space_for_lines(numLinesAlongWire, numTransverse, numLines, numPointsPerTurn, xArray, yArray, zArray, resolution); 



% Coordinates for plotting 
[xPlot, yPlot, zPlot] = extend_arrays_for_plots(numPoints, numThermalSubdivisions, ... 
    xArray, yArray, zArray, x1Array, x2Array, y1Array, y2Array, z1Array, z2Array, wCond, thCond, RSolEff); 

% Distance between the turns
xDiffArray = x2Array - x1Array;
yDiffArray = y2Array - y1Array;
zDiffArray = z2Array - z1Array; 

lenArray = sqrt(xDiffArray.^2 + yDiffArray.^2 + zDiffArray.^2); 
lenArrayAlongWire = lenArray(1:numLinesAlongWire); % Length of the conductor along turns, excluding transverse elements 
lenArrayAlongWireElements = lenArrayAlongWire*ones(1, numThermalSubdivisions)./numThermalSubdivisions; 

[MArray, LTotal] = generate_inductance_matrix(numLines, numLinesAlongWire, wCond, thCond, x1Array, x2Array, ... 
    y1Array, y2Array, z1Array, z2Array, mutualInductanceSpaceRatio, mutualInductanceSpaceRatioLimit, ignoreMutualInductancesTransverseElements);  

% Create matrix that keeps track of equations for lines and points
MLPInv = create_MLPInv_matrix(MArray, numLines, numPoints, pointIndex1, pointIndex2);

% ***** Biot-Savart calculations ***** 
% The special case of a uniform constant current I (current is out from the integral)

IArray = zeros(numLines, 1); % All lines (including transverse lines) can carry current
IArray(1:numLinesAlongWire) = I0; % Current is only induced along the line 
% IArray(1:3) = I0;
% IArray(end) = I0;

[BiotSavartMatrixX, BiotSavartMatrixY, BiotSavartMatrixZ, BLocalXArray, BLocalYArray, BLocalZArray] ... 
    = calc_biot_savart(numLines, x1Array, x2Array, y1Array, y2Array, z1Array, z2Array, ... 
    IArray, mutualInductanceSpaceRatio); 


% ***** Initial conditions of the magnet ***** 
temperatureArray = ones(numPoints - 1, numThermalSubdivisions)*initial_temperature;  
temperatureArray(floor(numPoints/2 + 0.5), floor(numThermalSubdivisions/2 + 1.5)) = 22; % Hotspot 22 K


% ***** Turn-to-turn resistance ***** 


arclengthCylinder = 2*pi*RSol;
widthWall = 0.005;          %m        


% ***** Helium gas cooling characteristics ***** 

A_helium = (wCond-heightSlot).*lenArrayAlongWireElements; %characteristic length for natural convection
L_helium = 0.5.*((wCond-heightSlot)+lenArrayAlongWireElements);
Nu = 4.36; %Nusselt number for steady laminar pipe flow


volumeArray = (lenArray(1:numLinesAlongWire)./numThermalSubdivisions ... 
    .*thCond*wCond)*ones(1, numThermalSubdivisions);
Acylinder = 2*2*pi*RSol*(wCond-heightSlot)*numWindings;
Aslots = 2*numTransverse*thSlot*thWallatShort;

Atotal = Acylinder+Aslots;
cooling_area_correction_factor = Atotal/Acylinder;



% The vector comprises all dI/dt followed by all voltage points, excluding
% point 1, which is per definition always at 0 V
V_local_element = zeros(numLinesAlongWire,1);
VGroundVector = zeros(numLinesAlongWire,1);

% ***** 

IExt = I0; % External current 

if writedata ==1
fprintf(fileID,'Input Parameters RePISoSi code: \n \n');
fprintf(fileID,'Initial Current: %.0f [A]\n', I0);
fprintf(fileID,'Radius solenoid: %.3f [m]\n', RSol);
fprintf(fileID,'Number of turns: %.1g [-]\n', numWindings);
fprintf(fileID,'Number of shorts per turn: %.2g [-]\n', N_shorts);
fprintf(fileID,'Number of line-elements between axial shorts: %.1g [-]\n', resolution);
fprintf(fileID,'Thickness of stabiliser wall: %.2g [mm]\n', 1e3*thCond);
fprintf(fileID,append('Tape type:',HTStapeName,'\n'));
fprintf(fileID,'Number of tapes [-]: %.1g\n', N_tapes);
fprintf(fileID,'Initial temperature [K]: %.1f\n', initial_temperature);
fprintf(fileID,'Power of Heater [W]: %.1f\n', PHeater);
fprintf(fileID,'Heating time [s]: %.1f\n', theater);
fprintf(fileID,'deltaTmax (used for determination time step) [K]: %.1f\n', dTMax);
fprintf(fileID,'temperatureSubSteps [-]: %.1f\n', temperatureSubSteps);
fprintf(fileID,'Status Helium convection cooling (1 = on, 0 = off): %.1g\n', Helium_cooling);
fprintf(fileID,'Status Copper heat conduction cooling (1 = on, 0 = off): %.1g\n', copper_heat);
fprintf(fileID,'Status Copper conductivitiy (1 = on, 0 = off): %.1g\n \n \n', copper_conduction);

fprintf(fileID,'time [s], centerB [T], Pheater [W], TringTop [K], TringBottom [K], Efield [uV/m], Tmax [K] \n');

end
disp('Transient calculation'); 

maxIteration = 50000; % Number of iterations 
updateCounter = floor(maxIteration/100);

drawFigureAtIteration = 20; % How many iterations go by before a plot is made.

s.MLPInv = MLPInv;     

temperatureCalculationStartTime = now*24*3600;
temperatureCalculationFinishTime = now*24*3600;
materialsPropsCalculationStartTime = now*24*3600;
materialsPropsCalculationFinishTime = now*24*3600;
ODECalculationStartTime = now*24*3600;
ODECalculationFinishTime = now*24*3600;

plottingStartTime = now*24*3600;
plottingFinishTime = now*24*3600;


overallFinishTime = now*24*3600;
overallStartTime = now*24*3600;
overallFirstStartTime = now*24*3600;




centerB = drawFieldMap(RSol, len, numPoints, numLines, numLinesAlongWire,numTransverse, numThermalSubdivisions, xPlot, yPlot, zPlot, ... 
    x1Array, x2Array, y1Array, y2Array, z1Array, z2Array, BLocalXArray, BLocalYArray, BLocalZArray, ... 
    temperatureArray, mutualInductanceSpaceRatio, IArray); 


t = 0; % Initial time (s) 
        temperatureRingTop = temperatureRing;
        temperatureRingBottom = temperatureRing;
% Linearized model to calculate quench propagation 
for iterationIndex = 1:maxIteration       

    overallStartTime = now*24*3600;
	materialsPropsCalculationStartTime = now*24*3600;
    

        
    externaldIdtVector = zeros(numPoints - 1, 1);
    
    % Current flows out of the last point
    s.externalVoltagedIdtVector = zeros(numPoints - 1, 1);

    % Material property calculation		    
    
    %thermalConductivityArray = thermal_conductivity_Al10SiMg(temperatureArray); % (W/(m*K)) 
    %thermalConductivityArray = thermal_conductivity_al_alloy_5083(temperatureArray); % (W/(m*K)) 
  
    %resistivityNormalArray = calc_normal_Al10SiMg_resistivity(temperatureArray); % (ohm*m)
    %resistivityNormalArray = calc_normal_Al_alloy_resistivity(temperatureArray); % (ohm*m)

    rounded(:,:) = round(10*temperatureArray(:,:));
    kmatrix = zeros(size(temperatureArray));
    kmatrix = ktable(rounded(:,:));
    thermalConductivityArray = kmatrix;
    rhomatrix = zeros(size(temperatureArray));
    rhomatrix = rhotable(rounded(:,:));
    resistivityNormalArray = rhomatrix;


    % Materials calculation starts here
	
	

    heatCapacityArray= heat_capacity_al_alloy_5083(temperatureArray, volumeArray);




    BXYZarray = [BLocalXArray BLocalYArray BLocalZArray]; %matrix of magnetic field direction at each point every row with columns 1,2,3 describing x,y,z coordinates respectively
    XYZarray =  [xArray yArray zArray];                   %matrix of cartesian coordinates of each point in every row with column  1,2,3 describing x,y,z coordinates respectively
    %BArray = sqrt(BLocalXArray.^2+BLocalYArray.^2+BLocalZArray.^2);
    BArray = vecnorm(BXYZarray,2,2); %length of vector corresponding to magnitude B-field, same as above

    
    CaxistapeArray = [xArray yArray zeros(length(xArray),1)];

    theta_angleArray = acos(dot(BXYZarray(1:length(CaxistapeArray),:),CaxistapeArray,2)./(BArray(1:length(CaxistapeArray)).*vecnorm(CaxistapeArray,2,2))); %angle in rad between coordinate vector at each point on cylinder, but with z=0, so parallel to the c-axis of the rebco tape and the magnetic field 

    maxTemparray = max(temperatureArray,[],2); %always get the highest temperature in the finer temperature scale
    IcArray = N_tapes*tapewidth*1e3*SCthickness*1e3*parametrisation_fujikura(BArray(1:length(temperatureArray)),maxTemparray,theta_angleArray(1:length(temperatureArray)));
 %   IcArray(1:round(length(IcArray)/5)) = 15000; 
  %  IcArray(round(4*length(IcArray)/5):end) = 15000;
    IcArrayhistory(:,iterationIndex) = IcArray;
	
	% Materials props calculation stops here and temperature part starts
	materialsPropsCalculationFinishTime = now*24*3600;
	
	
	%disp(sprintf("Materials calculation time: %0.3f s", materialsPropsCalculationFinishTime-materialsPropsCalculationStartTime));   
	
	temperatureCalculationStartTime = now*24*3600;




    temperatureBottomArray =  temperatureArray(1:resolution:(numTransverse-1)*resolution,1); %temperature at bottom element of axial short
    temperatureTopArray = temperatureArray(numPointsPerTurn+1:resolution:(numTransverse-1)*resolution+numPointsPerTurn,1);
    temperatureshort = (temperatureBottomArray+temperatureTopArray)/2;
    axialshortResistance = heightSlot*rhotable(round(temperatureshort*10))/(thSlot*widthWall); % Turn-to-turn resistance (ohm)|8.267091125040324e-08 : used to be 1e-6
    axialshortResistance = transpose(axialshortResistance);

    RArrayTrans = axialshortResistance;
    RArrayTrans(end+1) = RArrayTrans(end);

    RNormalArray = resistivityNormalArray.*lenArrayAlongWire(1:numLinesAlongWire)./(ACond.*numThermalSubdivisions);			
    if copper_conduction == 1
        RNormalArray(1:round(length(RNormalArray)/5),:) =         ones(round(length(RNormalArray)/5),3).*rhoCu_nist(temperatureRingBottom.*ones(round(length(RNormalArray)/5),1),BArray(1:round(length(RNormalArray)/5)),RingRRR,1).*lenArrayAlongWire(1:numLinesAlongWire/5)./(crossSectionRing.*numThermalSubdivisions); %consider copper for stabilisation of first turn
        RNormalArray(round(1+4*length(RNormalArray)/5):end,:) =     ones(round(length(RNormalArray)/5),3).*rhoCu_nist(temperatureRingBottom.*ones(round(length(RNormalArray)/5),1), ...
        BArray(round(1+4*length(RNormalArray)/5):round(length(RNormalArray))),RingRRR,1).*lenArrayAlongWire(1:numLinesAlongWire/5)./(crossSectionRing.*numThermalSubdivisions); %consider copper for stabilisation of first turn
    end

    IAbsArray = abs(IArray(1:numLinesAlongWire))*ones(1, numThermalSubdivisions); % Total current (A)
    %E0 = 100*1E-6; %[V/m] critical E-field convention
    %u0_array = E0.*lenArrayAlongWire(1:numLinesAlongWire); %[V] critical voltage over element following from Efield convention
    %n_value_array = n_value(maxTemparray,HTStapeName);

    
  
% 

%V_local_element = MArray(1:numLinesAlongWire,1:numLinesAlongWire) * dIdt;
% for n_index = 1:numLinesAlongWire
%     %fun = @(x) abs(L_self(n_index).*dIdt(n_index) + u0_array(n_index).*(x./IcArray(n_index)).^n_value_array(n_index) - (IAbsArray(n_index) - x).*RNormalArray(n_index));
%     %fun = @(x) abs(u0_array(n_index).*(x./IcArray(n_index)).^n_value_array(n_index) - (IArray(n_index) - x).*RNormalArray(n_index) - (IArraytrans2(n_index) - x).*RArraytrans2(n_index));
%     fun = @(x) abs(u0_array(n_index).*(x./IcArray(n_index)).^n_value_array(n_index) - abs((IArray(n_index) - x).*RNormalArray(n_index)) - abs(VGroundVector(n_index)));
% 
%     x = fminsearch(fun,4000);
% ISCArray(n_index,1) = x; %amount of current in superconductor for equal voltage [A]
% end

    INormalArray = IAbsArray - IcArray;               % Critical state model
       % INormalArray = IAbsArray - ISCArray;               % Powerlaw
    INormalArray = INormalArray.*(INormalArray > 0);  % Current element (I)  
    VElementArray = INormalArray.*RNormalArray;       % Voltage element (V)
    PElementArray = VElementArray.*IAbsArray;         % Heating power element (W = J/s)




    IArrayTrans = IArray(numLinesAlongWire+1:end-1);
    


    
    PElementArrayTransverse = IArrayTrans.^2.*axialshortResistance;


    PElementArray(1:resolution:(numTransverse-1)*resolution,1) = PElementArray(1:resolution:(numTransverse-1)*resolution,1) + 0.5*PElementArrayTransverse; %half of heating to top
    PElementArray(numPointsPerTurn+1:resolution:(numTransverse-1)*resolution+numPointsPerTurn,1) = PElementArray(numPointsPerTurn+1:resolution:(numTransverse-1)*resolution+numPointsPerTurn,1) + 0.5*PElementArrayTransverse; %half of heating to bottom


    halfindex = round(length(PElementArray)/2);

    % if  t < 1.6 %heating 2s
    %    PElementArray(halfindex,:) = PElementArray(halfindex,:) + PHeater/3; %[W] Force heating
    % end

    Pheatersave = 0;
    if toffset(1) < t && t < toffset(1) + theater(1) %heating 2s
       PElementArray(halfindex,:) = PElementArray(halfindex,:) + PHeater/3; %[W] Force heating
       Pheatersave = PHeater;
    end



    toffsettemp = toffset(1) + theater(1) + toffset(2);

    if toffsettemp < t && t < toffsettemp + theater(2) %heating 5s
       PElementArray(halfindex,:) = PElementArray(halfindex,:) + PHeater/3; %[W] Force heating
       Pheatersave = PHeater;
    end


    toffsettemp = toffsettemp + theater(2) + toffset(3);

    if toffsettemp < t && t < toffsettemp + theater(3) %heater 10s
       PElementArray(halfindex,:) = PElementArray(halfindex,:) + PHeater/3; %[W] %heating 10s
       Pheatersave = PHeater;
    end

    toffsettemp = toffsettemp + theater(3) + toffset(4);


    if toffsettemp < t && t < toffsettemp + theater(4) %heating 15s
       PElementArray(halfindex,:) = PElementArray(halfindex,:) + PHeater/3; %[W] %heating 10s
       Pheatersave = PHeater;
    end
    toffsettemp = toffsettemp + theater(4) + toffset(5);

    if toffsettemp < t && t < toffsettemp + theater(5) %heating 20s
       PElementArray(halfindex,:) = PElementArray(halfindex,:) + PHeater/3; %[W] %heating 10s
       Pheatersave = PHeater;
    end



    if t> 600
        break
    end


    
    

    dtInt = 0;
    

    % ***** Longitudinal heat flow, inside current line *****
    for temperatureIterationIndex = 1:temperatureSubSteps		

        temperatureArrayRight = [temperatureArray(:, 2:numThermalSubdivisions) temperatureArray(:, numThermalSubdivisions)]; %right boundary of cell
        temperatureArrayRight(1:numPoints - 2, numThermalSubdivisions) = temperatureArray(2:numPoints-1, 1);				%now every element sees the element to the right	and the element point sees itself. There are more points than lines, so thats why there is -2 instead of -1 

        thermalConductivityArrayRight = [thermalConductivityArray(:, 2:numThermalSubdivisions) thermalConductivityArray(:, numThermalSubdivisions)];
        thermalConductivityArrayRight(1:numPoints - 2, numThermalSubdivisions) = thermalConductivityArray(2:numPoints-1, 1);		

        temperatureArrayLeft = [temperatureArray(:, 1) temperatureArray(:, 1:numThermalSubdivisions - 1)];
        temperatureArrayLeft(2:numPoints - 1, 1) = temperatureArray(1:numPoints - 2 ,numThermalSubdivisions);

        thermalConductivityArrayLeft = [thermalConductivityArray(:, 1) thermalConductivityArray(:, 1:numThermalSubdivisions - 1)];
        thermalConductivityArrayLeft(2:numPoints-1, 1) = thermalConductivityArray(1:numPoints - 2, numThermalSubdivisions);



        % temperatureBottomArray =  temperatureArray(1:resolution:numTransverse*resolution,round(size(temperatureArray,2)/2)); %temperature at bottom element of axial short
        % temperatureTopArray = temperatureArray(numPointsPerTurn:resolution:(numTransverse-1)*resolution+numPointsPerTurn,round(size(temperatureArray,2)/2));
        % 
        % thermalConductivityBottomArray = thermalConductivityArray(1:resolution:numTransverse*resolution,round(size(temperatureArray,2)/2)); %thermal conductivity at bottom element of axial short. Take the middle 
        % thermalConductivityTopArray = thermalConductivityArray(numPointsPerTurn:resolution:(numTransverse-1)*resolution+numPointsPerTurn,round(size(temperatureArray,2)/2));

        temperatureBottomArray =  temperatureArray(1:resolution:(numTransverse-1)*resolution,1); %temperature at bottom element of axial short
        temperatureTopArray = temperatureArray(numPointsPerTurn+1:resolution:(numTransverse-1)*resolution+numPointsPerTurn,1);

        thermalConductivityBottomArray = thermalConductivityArray(1:resolution:(numTransverse-1)*resolution,1); %thermal conductivity at bottom element of axial short. Take the middle 
        thermalConductivityTopArray = thermalConductivityArray(numPointsPerTurn+1:resolution:(numTransverse-1)*resolution+numPointsPerTurn,1);

% ***** Heat Flux ********
        QElementArray_temp = (temperatureArrayRight - temperatureArray).*0.5 ... 
            .*(thermalConductivityArray + thermalConductivityArrayRight).*ACond./lenArrayAlongWireElements; %each element sees the element to the right of it and the boundary sees itself. 
        
        QElementArray_temp = QElementArray_temp + (temperatureArrayLeft - temperatureArray).*0.5 ... 
            .*(thermalConductivityArray + thermalConductivityArrayLeft).*ACond./lenArrayAlongWireElements; %each element sees the element to the left of it and the boundary sees itself. 



        Qaxialarray = (temperatureTopArray - temperatureBottomArray) .* 0.5.* (thermalConductivityBottomArray + thermalConductivityTopArray).*thSlot.*thWallatShort./heightSlot;

        QringBottomarray = (temperatureArray(1:numPointsPerTurn,:) - temperatureRingBottom).*(thermalConductivityArray(1:numPointsPerTurn,:)).*lenArrayAlongWireElements(1:numPointsPerTurn,:)*thCond./(0.5.*wCond); %approximation that local thermal conductivity of aluminium is limiting the heat flow to copper and that copper stays at constant temperature. Also, the thermal conductivity is calculated as the conductivity as if the heat would travel axially through an aluminium thickness of half a turn. Contact resistance to copper is neglected.

        QringToparray = (temperatureArray(end+1-numPointsPerTurn:end,:) - temperatureRingTop).*(thermalConductivityArray(end+1-numPointsPerTurn:end,:)).*lenArrayAlongWireElements(end+1-numPointsPerTurn:end,:)*thCond./(0.5.*wCond);

        Qheliumcoolingarray(:,:) = 0;
        if Helium_cooling == 1
            k_heliumarray = k_helium*ones(size(temperatureArray,1),size(temperatureArray,2)); %helium temperature constant

            h_helium_array = Nu*k_heliumarray./L_helium; %solving heat transfer coefficient from Nusselt number

            R_thermal = 0.5*thCond./(thermalConductivityArray)+1./h_helium_array; %Thermal resistance, including heat conductivity

            Qheliumcoolingarray = 2.*cooling_area_correction_factor*A_helium.*(temperatureArray-temperature_helium)./R_thermal; %(area from 2 sides)




        end


        % 
        % QElementArray_temp(1:resolution:numTransverse*resolution,round(size(temperatureArray,2)/2)) = QElementArray_temp(1:resolution:numTransverse*resolution,round(size(temperatureArray,2)/2)) + Qaxialarray; %add heat flow for bottom part sign convention 
        % QElementArray_temp(numPointsPerTurn:resolution:(numTransverse-1)*resolution+numPointsPerTurn,round(size(temperatureArray,2)/2)) = QElementArray_temp(numPointsPerTurn:resolution:(numTransverse-1)*resolution+numPointsPerTurn,round(size(temperatureArray,2)/2)) - Qaxialarray; %add heat flow for bottom part sign convention 

        
        QElementArray_temp(1:resolution:(numTransverse-1)*resolution,1) = QElementArray_temp(1:resolution:(numTransverse-1)*resolution,1) + Qaxialarray; %add heat flow for bottom part sign convention 
        QElementArray_temp(numPointsPerTurn+1:resolution:(numTransverse-1)*resolution+numPointsPerTurn,1) = QElementArray_temp(numPointsPerTurn+1:resolution:(numTransverse-1)*resolution+numPointsPerTurn,1) - Qaxialarray; %add heat flow for bottom part sign convention 

        

        if copper_heat == 1
%%%Add copper at the bottom
        QElementArray_temp(1:numPointsPerTurn,:) = QElementArray_temp(1:numPointsPerTurn,:) -QringBottomarray; % - since QringBottomarray is + if temperature of ring smaller than temperature of cylinder

%%%Add copper at the top
        QElementArray_temp(end+1-numPointsPerTurn:end,:) = QElementArray_temp(end+1-numPointsPerTurn:end,:) - QringToparray; % - since QringToparray is + if temperature of ring smaller than temperature of cylinder

        end


        QElementArray = QElementArray_temp - Qheliumcoolingarray;

        %dTdtElementArray = (PElementArray + QElementArray)./(heatCapacityArray.*volumeArray); Wrong IMO, since       %heat capacity already includes volume
        dTdtElementArray = (PElementArray + QElementArray)./(heatCapacityArray); 
        dtLocal = dTMax./max(max(abs(dTdtElementArray)))/temperatureSubSteps;
        

        heatCapacityCoppertop = CpCu_nist(temperatureRingTop,massRing);
         heatCapacityCopperbottom = CpCu_nist(temperatureRingBottom,massRing);

        dTdtCoppertop = sum(sum(QringToparray))/heatCapacityCoppertop;
        dTdtCopperBottom = sum(sum(QringBottomarray))/heatCapacityCopperbottom;


        if max(max(abs(dTdtElementArray)))/temperatureSubSteps == 0
            disp('divide by 0')
        end

        %nancheck(iterationIndex,temperatureIterationIndex) = dtLocal ;%bugfixing
        dtInt = dtInt + dtLocal;
        %dtIntarray(iterationIndex,temperatureIterationIndex) = dtInt; %bugfixing
        temperatureArray = temperatureArray + dtLocal*dTdtElementArray;
        temperatureRingBottom = temperatureRingBottom +dtLocal*dTdtCopperBottom;
        temperatureRingTop = temperatureRingTop + dtLocal*dTdtCoppertop;

    end
    temperatureArraystored(:,:,iterationIndex) = temperatureArray;
    dt = dtInt;

    temperatureRingTophistory(:,iterationIndex) =   temperatureRingTop;
    temperatureRingBottomhistory(:,iterationIndex) = temperatureRingBottom;
    temperaturemax(iterationIndex) = max(temperatureArray(:));
    
	
	% Temperature part stops here and ODE part starts
	temperatureCalculationFinishTime = now*24*3600;
	%disp(sprintf("Temperature calculation time: %0.3f s", temperatureCalculationFinishTime-temperatureCalculationStartTime));   
	ODECalculationStartTime = now*24*3600;
	
	

    % Construction of vector b in Av = b
    % The vector comprises all dI/dt followed by all voltage points, excluding
    % point 1, which is per definition always at 0 V. 
    
    s.RArrayTrans = RArrayTrans;
    s.IcArray = IcArray;   %critical current of superconductor
    %s.ISCArray = ISCArray;  %current in superocnductor
    s.externalVoltagedIdtVector = s.externalVoltagedIdtVector;
    s.bVector = [];
    s.VGroundVector = [];
    s.numLines = numLines;
    s.numLinesAlongWire = numLinesAlongWire;
    s.numTransverse = numTransverse;
    s.numPoints = numPoints;
    s.lenArrayAlongWire = lenArrayAlongWire;
    s.signArray = zeros(numLinesAlongWire, 1);
    s.INormalArray = zeros(numLinesAlongWire, 1);
    s.RNormalArray = RNormalArray;
    s.numThermalSubdivisions = numThermalSubdivisions;

    IArray0 = IArray;

    tspan = [t, (t + dt)];      
    
    myfun = @(t, IArray)myODE(t, IArray, s); 
    % [tArchive, IArrayArchive] = ode23tb(myfun, tspan, IArray0); % Function ode23tb solves stiff differential equations, low order method 
 %   [tArchive, IArrayArchive] = ode23s(myfun, tspan, IArray0); %GOOD ONE
        [tArchive, IArrayArchive] = ode15s(myfun, tspan, IArray0); %GOOD ONE

    %[tArchive, IArrayArchive] = ode89(myfun, tspan, IArray0); 
    t = tArchive(length(tArchive));
    IArray = IArrayArchive(length(tArchive), :)';        


    tArrayhistory(iterationIndex) = t;
    IArrayhistory(:,iterationIndex) = IArray;
    VGroundVector = myODEVoltages(t, IArray, s);
    % V_local_element = diff(VGroundVector);
    % V_local_element(end+1) = 0;


    Qheliumcoolingarrayhistory(:,:,iterationIndex) = Qheliumcoolingarray;
    VGroundVectorhistory(:,iterationIndex) = VGroundVector;
    Efield(iterationIndex) = 1e6*(VGroundVector(round(length(VGroundVector)*1/10))-VGroundVector(round(length(VGroundVector)*9/10)))./(4*2*pi*0.12); % Be very careful with these variables

    % ODE part finishes here and plotting part starts
	ODECalculationFinishTime = now*24*3600;	
	%disp(sprintf("ODE calculation time: %0.3f s", ODECalculationFinishTime-ODECalculationStartTime));   
	
	
	plottingStartTime = now*24*3600;

    % Local magnetic flux density (T)
    BLocalXArray = BiotSavartMatrixX*IArray; 
    BLocalYArray = BiotSavartMatrixY*IArray; 
    BLocalZArray = BiotSavartMatrixZ*IArray; 

    if(mod(iterationIndex+1,drawFigureAtIteration) == 0)
        
        centerB = drawFieldMap(RSol, len, numPoints, numLines, numLinesAlongWire, numTransverse,numThermalSubdivisions, xPlot, yPlot, zPlot, ... 
            x1Array, x2Array, y1Array, y2Array, z1Array, z2Array, BLocalXArray, BLocalYArray, BLocalZArray, ... 
            temperatureArray, mutualInductanceSpaceRatio, IArray); 
        title(sprintf('Iteration = %0.0f, t = %0.3g s, I_{ext} = %0.2f A', iterationIndex, t, IExt));
        pause(0.1) 
        figure(1);
        hold on

        if writevideo ==1
        frame = getframe(gcf);
        writeVideo(v,frame);
        end

        hold off
            centerBhistory(iterationIndex) = centerB;
    figure(2)
    plot(tArrayhistory(centerBhistory>0),centerBhistory(centerBhistory>0),Color='b',Marker='.',MarkerSize=10,LineStyle='-')
    xlabel('t [s]')
    ylabel('B_{center} [T]')
    grid on

    figure(3)
    plot(tArrayhistory,Efield,Color='b',Marker='.',MarkerSize=10,LineStyle='-')
    xlabel('t [s]')
    ylabel('E_{field} [\muV/m]')
    grid on
    end 
    
    

    %dIdt = myODE(t, IArray, s);
	
	% Plotting part finishes here
	plottingFinishTime = now*24*3600;	
	%disp(sprintf("Plotting time: %0.3f s", plottingFinishTime-plottingStartTime));   
	overallFinishTime = now*24*3600;
	%disp(sprintf("Overall time: %0.3f s", overallFinishTime-overallStartTime));   
	
	
	    %Calculate fraction of time spent in each operation
    fracMatCalc = 100*(materialsPropsCalculationFinishTime - materialsPropsCalculationStartTime)./(overallFinishTime - overallStartTime); 
	fracTSteps = 100*(temperatureCalculationFinishTime - temperatureCalculationStartTime)./(overallFinishTime - overallStartTime);
	fracODE = 100*(ODECalculationFinishTime - ODECalculationStartTime)./(overallFinishTime - overallStartTime); 
	fracPlotting = 100*(plottingFinishTime - plottingStartTime)./(overallFinishTime - overallStartTime);    
    

    disp(sprintf("Index: %0.0f, time expired: %0.1f s, tsim %0.6g s, last iteration %0.1f s (%0.3f%% in material calc, %0.3f%% in TSteps, %0.3f%% in ODE, %0.3f%% plotting)",iterationIndex, overallFinishTime-overallFirstStartTime, t, overallFinishTime-overallStartTime, fracMatCalc, fracTSteps, fracODE, fracPlotting));   

    
    if writedata ==1
    fprintf(fileID,'%0.6g, %0.6f, %0.3f, %0.3f, %0.3f, %0.6f,  %0.4f \n',[tArrayhistory(1,end) centerB Pheatersave temperatureRingTop temperatureRingBottom Efield(1,end) temperaturemax(1,end)]);
    end
	

end 

if writedata ==1
fclose(fileID);
end

if writevideo ==1
close(v);
end
save(append(Description,string(datetime('now','TimeZone','local','Format','d-MMM-y')),'_',string(datetime('now','TimeZone','local','Format','HH.mm.ss')),'.mat'))
toc       
    