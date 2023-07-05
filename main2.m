% ***** Quench analysis of partially insulated coil *****
% 
%
% Program written by M. Mentink and A. Vaskuri (Oct. 27, 2022). 
%
% *****


clear all; close all; clc; 

set(0, 'defaultTextFontName', 'times', 'defaultTextFontSize', 14);
set(0, 'defaultAxesFontName', 'times', 'defaultAxesFontSize', 14);

format long

% ***** Solenoid properties ***** 
RSol = 0.2;                 % Solenoid radius (m)

numPointsPerTurn = 30;       % Sets the resolution 
wCond = 0.004;              % Conductor width (m) 
thCond = 0.0038;             % Effective conductor thickness (m), 
                             % 3.8 mm = 1.5 mm * 5 mm + 4.5 mm * 3 mm + 1.5 mm * 5 mm)/(1.5 mm + 4.5 mm + 1.5 mm)   

thSlot = 0.001;              % Slot thickness (m)
ACond = wCond.*thCond;       % Cross-sectional area of the conductor (m)

numWindings = 20;            % Number of turn forming a solenoid (dimensionless)
len = numWindings*(wCond + thSlot);     % Length of a solenoid (m)

ignoreMutualInductancesTransverseElements = 1; % (yes = 1, no = 0)
numThermalSubdivisions = 5;  % Number of subdivisions 
temperatureSubSteps = 100;   % Number of temperature increments 

dTMax = 1;                   % Maximum temperature difference (K)

I0 = 1000;                    % Driving current (A) 
% ***** 


% ***** Superconducting cable properties ***** 
Tc = 77;   % Critical temperature (K) 
Ic0 = 2000; % Critical current near 0 K (A) 
% ***** 


% ***** Look up tables ***** 
Tinitialise = 0.1:0.1:400;
L = 2.44E-8; %[V^2K^-2]
ktable = 1E6.*Tinitialise.*L./BlochGrun(Tinitialise,0.2594,484,0.026,5);
rhotable = BlochGrun(Tinitialise,0.2594,484,0.026,5)/1E6;
% ***** 

%%

% ***** Initialize mesh for FEM ***** 
pointAngle = 2*pi/numWindings;  % Angle between the points (rad)
RSolMin = RSol*cos(pointAngle); % Minimum distance from the center due to finite number of points (m)
RSolAvg = 0.5*(RSolMin + RSol); % Average distance (m)

RSolEff = RSol^2/RSolAvg;       % Effective distance (m)
mutualInductanceSpaceRatio = 10; 
mutualInductanceSpaceRatioLimit = 100;

numPoints = floor(numPointsPerTurn*numWindings + 1);            % Total number of 3D points
numThermalElements = (numPoints - 1)*numThermalSubdivisions;    % Total number of thermal elements

% Generate line elements 
[xArray, yArray, zArray] = ... 
    generate_line_elements(numPoints, numPointsPerTurn, RSolEff, numWindings, wCond); 

numLinesAlongWire = numPoints - 1; 
numTransverse = floor((numWindings - 1)*numPointsPerTurn + 1); 
numTransverse = max(numTransverse, 0); 
numLines = numLinesAlongWire + numTransverse; 

% Gather space for lines
[x1Array, x2Array, y1Array, y2Array, z1Array, z2Array, pointIndex1, pointIndex2] = ... 
    gather_space_for_lines(numLinesAlongWire, numTransverse, numLines, numPointsPerTurn, xArray, yArray, zArray); 

% Coordinates for plotting 
[xPlot, yPlot, zPlot] = extend_arrays_for_plots(numPoints, numThermalSubdivisions, ... 
    xArray, yArray, zArray, x1Array, x2Array, y1Array, y2Array, z1Array, z2Array, wCond, thCond, RSolEff); 

% Distance between the turns             
xDiffArray = x2Array - x1Array;
yDiffArray = y2Array - y1Array;
zDiffArray = z2Array - z1Array; 

lenArray = sqrt(xDiffArray.^2 + yDiffArray.^2 + zDiffArray.^2); 
lenArrayAlongWire = lenArray(1:numLinesAlongWire); % Length of the conductor per turn 
lenArrayAlongWireElements = lenArrayAlongWire*ones(1, numThermalSubdivisions)./numThermalSubdivisions; 

[MArray, LTotal] = generate_inductance_matrix(numLines, numLinesAlongWire, wCond, thCond, x1Array, x2Array, ... 
    y1Array, y2Array, z1Array, z2Array, mutualInductanceSpaceRatio, mutualInductanceSpaceRatioLimit, ignoreMutualInductancesTransverseElements);  

% Create matrix that keeps track of equations for lines and points
MLPInv = create_MLPInv_matrix(MArray, numLines, numPoints, pointIndex1, pointIndex2);

% ***** Biot-Savart calculations ***** 
% The special case of a uniform constant current I (current is out from the integral)

IArray = zeros(numLines, 1); % All lines (including transverse lines) can carry current
IArray(1:numLinesAlongWire) = I0; % Current is only induced along the line 

[BiotSavartMatrixX, BiotSavartMatrixY, BiotSavartMatrixZ, BLocalXArray, BLocalYArray, BLocalZArray] ... 
    = calc_biot_savart(numLines, x1Array, x2Array, y1Array, y2Array, z1Array, z2Array, ... 
    IArray, mutualInductanceSpaceRatio); 


% ***** Initial conditions of the magnet ***** 
initial_temperature = 4.2; % Initial temperature of 4 K
%initial_temperature = 77; % Initial temperature of 77 K
temperatureArray = ones(numPoints - 1, numThermalSubdivisions)*initial_temperature;  
%temperatureArray(floor(numPoints/2 + 0.5), floor(numThermalSubdivisions/2 + 1.5)) = 10; % Hotspot 10 K
temperatureArray(floor(numPoints/2 + 0.5), floor(numThermalSubdivisions/2 + 1.5)) = 10; % Hotspot 100 K


% ***** Turn-to-turn resistance ***** 
% The vector comprises all dI/dt followed by all voltage points, excluding
% point 1, which is per definition always at 0 V. 
turnToTurnResistance = 1e-6; % Turn-to-turn resistance (ohm)

RArrayTrans = ones(numTransverse, 1)*turnToTurnResistance/numPointsPerTurn; % Transverse resistance array (ohm) 
%RArrayTrans = turnToTurnResistance/numPointsPerTurn; % Transverse resistance array (ohm) 
%%

volumeArray = (lenArray(1:numLinesAlongWire)./numThermalSubdivisions ... 
    .*thCond*wCond)*ones(1, numThermalSubdivisions);

% ***** 

IExt = I0; % External current 
%IArray(1:numLinesAlongWire) = I0; % Initial current 

disp('Transient calculation'); 

maxIteration = 300; % Number of iterations 
updateCounter = floor(maxIteration/100);

drawFigureAtIteration = 1; % (yes 1, no 0)

s.MLPInv = MLPInv;      

timeOdeNext = now*24*3600;
timeOdePrevious = timeOdeNext;
timeTemperatureNext = now*24*3600;
timeTemperaturePrevious = timeTemperatureNext;
timeMaterialCalculationNext = now*24*3600;
timeMaterialCalculationPrevious = timeMaterialCalculationNext;
timeLastIteration = now*24*3600;
timeNextIteration = now*24*3600;

pause(0.1); 

drawFieldMap(RSol, len, numPoints, numLines, numLinesAlongWire, numThermalSubdivisions, xPlot, yPlot, zPlot, ... 
    x1Array, x2Array, y1Array, y2Array, z1Array, z2Array, BLocalXArray, BLocalYArray, BLocalZArray, ... 
    temperatureArray, mutualInductanceSpaceRatio, IArray); 

t = 0; % Initial time (s) 

% Linearized model to calculate quench propagation 
for iterationIndex = 1:maxIteration       
    timeNextIteration = now*24*3600; % Time stamp in seconds 

    fracODE = 100*(timeOdeNext - timeOdePrevious)./(timeNextIteration - timeLastIteration);
    fracTSteps = 100*(timeTemperatureNext - timeTemperaturePrevious)./(timeNextIteration - timeLastIteration);
    fracMatCalc = 100*(timeMaterialCalculationNext - timeMaterialCalculationPrevious)./(timeNextIteration - timeLastIteration); 

    disp(sprintf("Index: %0.0f, time expired: %0.1f s, tsim %0.6g s, last iteration %0.1f s (%0.0f%% in ODE, %0.0f%% in TSteps, %0.0f%% in material calc)",iterationIndex, toc, t, timeNextIteration-timeLastIteration, fracODE, fracTSteps, fracMatCalc));   
    timeLastIteration = timeNextIteration;
    
    externaldIdtVector = zeros(numPoints - 1, 1);
    
    % Current flows out of the last point
    s.externalVoltagedIdtVector = zeros(numPoints - 1, 1);

    % Material property calculation		
    timeMaterialCalculationPrevious = now*24*3600;  % (s) 
    
    %thermalConductivityArray = thermal_conductivity_Al10SiMg(temperatureArray); % (W/(m*K)) 
    thermalConductivityArray = thermal_conductivity_al_alloy_5083(temperatureArray); % (W/(m*K)) 
  
    %resistivityNormalArray = calc_normal_Al10SiMg_resistivity(temperatureArray); % (ohm*m)
    resistivityNormalArray = calc_normal_Al_alloy_resistivity(temperatureArray); % (ohm*m)

%     rounded(:,:) = round(10*temperatureArray(:,:));
%     kmatrix = zeros(size(temperatureArray));
%     kmatrix = ktable(rounded(:,:));
%     thermalConductivityArray = kmatrix;
%     rhomatrix = zeros(size(temperatureArray));
%     rhomatrix = rhotable(rounded(:,:));
%     resistivityNormalArray = rhomatrix;

    heatCapacityArray = heat_capacity_al_alloy_5083(temperatureArray, sum(volumeArray(:))); 
    checkresistivityarray(:,:,iterationIndex) = resistivityNormalArray;
    checkheatconductivityarray(:,:,iterationIndex) = thermalConductivityArray;

    IcArray = calc_HTS_critical_current(temperatureArray, Ic0, Tc); % (A) 

    turnToTurnResistance = turnToTurnResistance; % TODO: Must be updated if temperature changes  
    
    RArrayTrans = turnToTurnResistance/numPointsPerTurn; % Transverse resistance array (ohm) 
    
    RNormalArray = resistivityNormalArray.*lenArrayAlongWire./(ACond.*numThermalSubdivisions);			

    IAbsArray = abs(IArray(1:numLinesAlongWire))*ones(1, numThermalSubdivisions); % Total current (A)
    INormalArray = IAbsArray - IcArray;               % Critical state model
    INormalArray = INormalArray.*(INormalArray > 0);  % Current element (I)
    VElementArray = INormalArray.*RNormalArray;       % Voltage element (V)
    PElementArray = VElementArray.*IAbsArray;         % Heating power element (W = J/s)

    timeMaterialCalculationNext = now*24*3600;

    % Thermal part
    timeTemperaturePrevious = now*24*3600;	

    dtInt = 0;

    % ***** Longitudinal heat flow, inside current line *****
    for temperatureIterationIndex = 1:temperatureSubSteps		
        temperatureArrayRight = [temperatureArray(:, 2:numThermalSubdivisions) temperatureArray(:, numThermalSubdivisions)];
        temperatureArrayRight(1:numPoints - 2, numThermalSubdivisions) = temperatureArray(2:numPoints - 1, 1);					
        thermalConductivityArrayRight = [thermalConductivityArray(:, 2:numThermalSubdivisions) thermalConductivityArray(:, numThermalSubdivisions)];
        thermalConductivityArrayRight(1:numPoints - 2, numThermalSubdivisions) = thermalConductivityArray(2:numPoints - 1, 1);		

        temperatureArrayLeft = [temperatureArray(:, 1) temperatureArray(:, 1:numThermalSubdivisions - 1)];
        temperatureArrayLeft(2:numPoints - 1, 1) = temperatureArray(1:numPoints - 2 ,numThermalSubdivisions);
        thermalConductivityArrayLeft = [thermalConductivityArray(:, 1) thermalConductivityArray(:, 1:numThermalSubdivisions - 1)];
        thermalConductivityArrayLeft(2:numPoints-1, 1) = thermalConductivityArray(1:numPoints - 2, numThermalSubdivisions);

        QElementArray_temp = (temperatureArrayRight - temperatureArray).*0.5 ... 
            .*(thermalConductivityArray + thermalConductivityArrayRight).*ACond./lenArrayAlongWireElements; 
        
        QElementArray = QElementArray_temp + (temperatureArrayLeft - temperatureArray).*0.5 ... 
            .*(thermalConductivityArray + thermalConductivityArrayLeft).*ACond./lenArrayAlongWireElements; 

        dTdtElementArray = (PElementArray + QElementArray)./(heatCapacityArray.*volumeArray);

        dtLocal = dTMax./max(max(abs(dTdtElementArray)))/temperatureSubSteps;
        if max(max(abs(dTdtElementArray)))/temperatureSubSteps == 0
            disp('divide by 0')

        end

        nancheck(iterationIndex,temperatureIterationIndex) = dtLocal ;%bugfixing
        dtInt = dtInt + dtLocal;
        dtIntarray(iterationIndex,temperatureIterationIndex) = dtInt; %bugfixing
        temperatureArray = temperatureArray + dtLocal*dTdtElementArray;
        
    end
    temperatureArraystored(:,:,iterationIndex) = temperatureArray;
    dt = dtInt;
    % *****

    % Material property calculation
    timeTemperatureNext = now*24*3600; 

    % Construction of vector b in Av = b
    % The vector comprises all dI/dt followed by all voltage points, excluding
    % point 1, which is per definition always at 0 V. 
    
    s.RArrayTrans = RArrayTrans;
    s.IcArray = IcArray;      
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
    
    timeOdePrevious = timeOdeNext;
    myfun = @(t, IArray)myODE(t, IArray, s); 
    % [tArchive, IArrayArchive] = ode23tb(myfun, tspan, IArray0); % Function ode23tb solves stiff differential equations, low order method 
    [tArchive, IArrayArchive] = ode23s(myfun, tspan, IArray0); 
    t = tArchive(length(tArchive));
    IArray = IArrayArchive(length(tArchive), :)';        
    VGroundVector = myODEVoltages(t, IArray, s);

    timeOdeNext = now*24*3600;

    % Local magnetic flux density (T)
    BLocalXArray = BiotSavartMatrixX*IArray; 
    BLocalYArray = BiotSavartMatrixY*IArray; 
    BLocalZArray = BiotSavartMatrixZ*IArray; 

    if(mod(iterationIndex,drawFigureAtIteration) == 0)
        drawFieldMap(RSol, len, numPoints, numLines, numLinesAlongWire, numThermalSubdivisions, xPlot, yPlot, zPlot, ... 
            x1Array, x2Array, y1Array, y2Array, z1Array, z2Array, BLocalXArray, BLocalYArray, BLocalZArray, ... 
            temperatureArray, mutualInductanceSpaceRatio, IArray); 
        title(sprintf('Iteration = %0.0f, t = %0.3f s, I_{ext} = %0.2f A', iterationIndex, t, IExt));
        pause(0.1) 
    end 
end 

    
toc       
    