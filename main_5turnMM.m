% =========================================================================
%  main_5turnMM.m  --  RePISoSi: Quench analysis of a partially-insulated
%                                 ReBCO solenoid (canonical entry point).
%
%  Author : J.L. Van den Eijnden, 2026
%           Original implementation by M. Mentink (Oct. 2022),
%           extended by J.L. Van den Eijnden and A. Vaskuri  (Oct. 2023).
%  License: MIT (see LICENSE in the repository root)
%
%  Purpose
%  -------
%  Couples a 3-D Biot-Savart self/mutual-inductance network model of an
%  HTS solenoid with a finite-volume thermal solver and a stiff ODE for
%  the line-element currents (ode15s). Models partially-insulated turn-to-
%  turn shorts, helium gas convection, and copper current leads.
%
%  Outputs (toggled via writedata / writevideo at the top of the
%  CONFIGURATION section):
%    * Live 3-D field-map figure
%    * <Description>.mp4 video
%    * <Description><timestamp>.txt   per-iteration log
%    * <Description><timestamp>.mat   final workspace dump
%
%  Usage
%  -----
%      addpath(genpath(pwd));   % once per MATLAB session
%      run('main_5turnMM.m');
%
%  Edit the CONFIGURATION block below to change the run.
% =========================================================================

clear; close all; clc;

set(0, 'defaultTextFontName', 'times', 'defaultTextFontSize', 14);
set(0, 'defaultAxesFontName', 'times', 'defaultAxesFontSize', 14);
format long;

tic;

%% ===================== CONFIGURATION =====================================
% ---- Run identification & I/O toggles -----------------------------------
HTStapeName = 'Fujikura FESC-04 4 mm wide';   % Tape used (informational)
Description = '15turn4431A_80VHeater_allon';  % TODO: parameterize per run
writedata   = 1;                              % 1 = write CSV-style log file
writevideo  = 0;                              % 1 = write MP4 of field map

% ---- Physics flags ------------------------------------------------------
Helium_cooling     = 1;   % 1 = on, 0 = off
copper_heat        = 1;   % copper-lead heat conduction (1=on, 0=off)
copper_conduction  = 1;   % copper-lead electrical conduction (1=on, 0=off)

% ---- Driving current ----------------------------------------------------
I0 = 4431;        % Driving current (A)

% ---- Heater pulse profile (5 pulses) -----------------------------------
% toffset(k) = wait time before pulse k (s); theater(k) = duration of pulse k (s)
toffset = [8.071, 48.795, 70.456, 84.040, 69.285, 100];   % s
theater = [1.633,  5.478, 10.304, 15.420, 19.826];        % s
RHeater = 350;        % heater resistance (ohm)
VHeater = 80;         % heater voltage (V)
PHeater = VHeater.^2 / RHeater;   % heater power (W)

% ---- Time-stepping control ----------------------------------------------
dTMax               = 1;     % Max temperature change per substep (K)
temperatureSubSteps = 200;   % Number of thermal substeps per magnetic step
maxIteration        = 50000; % Number of magnetic-network iterations
drawFigureAtIteration = 20;  % Plot field map every N iterations
simulation_time_limit = 600; % Hard stop at this simulated time (s)

% ---- Initial temperatures ----------------------------------------------
initial_temperature = 20;                       % Coil initial temperature (K)
temperature_helium  = initial_temperature;      % He bath temperature (K)
k_helium            = conductivity_helium(temperature_helium);

% ---- Solenoid geometry --------------------------------------------------
RSol         = 0.12;       % Solenoid radius (m)
N_shorts     = 19;         % Number of axial shorts per turn
arclengthSlot = 4e-2;      % Arclength of slot (m, informational)
thSlot        = 3.7e-3;    % Arclength of short (m)

% Number of line elements between two axial shorts (must be even).
resolution = 4;
resolution = round(resolution/2) * 2;

% Off-by-one-turn axial-short placement maximises the axial resistance.
numPointsPerTurn = N_shorts*resolution - resolution/2;

wCond        = 0.0085;     % Conductor width (m)
thReBCO      = 2.2e-6;     % ReBCO layer thickness (m)
% Effective Al-stabiliser thickness (incl. gap):
%   3.8 mm = (1.5*5 + 4.5*3 + 1.5*5) / (1.5 + 4.5 + 1.5)
thCond       = 0.0038;     % Effective conductor thickness (m)
thWallatShort = 5e-3;      % Wall thickness between slots (m)
heightSlot   = 1e-3;       % Slot height (m)
ACond        = wCond .* thCond;        % Conductor cross-sectional area (m^2)

numWindings  = 5;                      % Number of turns
len          = numWindings * (wCond + thSlot);   % Solenoid axial length (m)

ignoreMutualInductancesTransverseElements = 0;
numThermalSubdivisions = 3;            % Subdivisions per element (thermal)

% ---- Current-lead (copper ring) properties ------------------------------
NRings      = 2;
DouterLead  = 0.307;        % Outer diameter (m)
DinnerLead  = 0.208;        % Inner diameter (m)
ARing       = pi*(DouterLead^2 - DinnerLead^2)/4;     % Annular area (m^2)
thRing      = 0.01;                                    % Ring thickness (m)
crossSectionRing = thRing*(DouterLead - DinnerLead);   % Radial cross-section (m^2)
rhoCopper   = 8960;                                    % Copper density (kg/m^3)
massRing    = ARing * thRing * rhoCopper;              % Ring mass (kg)
temperatureRing = initial_temperature;                 % Initial ring T (K)
% RRR ~110 (normal Cu), ~45 (non-annealed OFE), 160 used here for annealed OFE.
RingRRR     = 160;

% ---- Superconducting tape parameters ------------------------------------
N_tapes     = 4;
tapewidth   = 4e-3;     % m
SCthickness = 2.2e-6;   % m

% ---- Mutual-inductance numerical parameters -----------------------------
mutualInductanceSpaceRatio      = 10;
mutualInductanceSpaceRatioLimit = 1000;

%% ===================== MATERIAL LOOKUP TABLES ============================
% Lookup tables are indexed by round(10*T), so T resolution is 0.1 K.
Tinitialise = 0.1:0.1:400;
L_lorenz    = 2.44e-8;     % Lorenz number (V^2 K^-2)
% Bloch-Grueneisen fit parameters for Al-10SiMg taken from prior calibration.
ktable   = 1e6 .* Tinitialise .* L_lorenz ./ ...
           BlochGrun(Tinitialise, 0.24233, 432.43,    0.0247015, 5);
rhotable = BlochGrun(Tinitialise, 0.24233, 432.43717, 0.0247015, 5) / 1e6;

%% ===================== MESH GENERATION ===================================
% Effective radius accounts for finite point count (chord vs. arc).
pointAngle = 2*pi/numPointsPerTurn;
RSolMin    = RSol * cos(pointAngle);
RSolAvg    = 0.5 * (RSolMin + RSol);
RSolEff    = RSol^2 / RSolAvg;

numPoints           = floor(numPointsPerTurn*numWindings + 1);
numThermalElements  = (numPoints - 1) * numThermalSubdivisions;
numLinesAlongWire   = numPoints - 1;

% Number of axial shorts. The "-round(numWindings/2)" term accounts for the
% deliberate off-by-one-turn pattern: each round leaves ~0.5 turn unconnected.
numTransverse = floor((numWindings - 1)*N_shorts + 1) - round(numWindings/2);
numTransverse = max(numTransverse, 0);
numLines      = numLinesAlongWire + numTransverse;

dIdt        = zeros(numLinesAlongWire, 1);
IArrayTrans = zeros(numTransverse - 1, 1);

% --- Helix line elements --------------------------------------------------
[xArray, yArray, zArray] = generate_line_elements( ...
    numPoints, numPointsPerTurn, RSolEff, numWindings, wCond);

[x1Array, x2Array, y1Array, y2Array, z1Array, z2Array, ...
 pointIndex1, pointIndex2] = gather_space_for_lines( ...
    numLinesAlongWire, numTransverse, numLines, numPointsPerTurn, ...
    xArray, yArray, zArray, resolution);

% --- Plotting coordinates -------------------------------------------------
[xPlot, yPlot, zPlot] = extend_arrays_for_plots(numPoints, ...
    numThermalSubdivisions, xArray, yArray, zArray, ...
    x1Array, x2Array, y1Array, y2Array, z1Array, z2Array, ...
    wCond, thCond, RSolEff);

% --- Element lengths ------------------------------------------------------
xDiffArray = x2Array - x1Array;
yDiffArray = y2Array - y1Array;
zDiffArray = z2Array - z1Array;
lenArray   = sqrt(xDiffArray.^2 + yDiffArray.^2 + zDiffArray.^2);
lenArrayAlongWire         = lenArray(1:numLinesAlongWire);
lenArrayAlongWireElements = lenArrayAlongWire * ones(1, numThermalSubdivisions) ...
                            ./ numThermalSubdivisions;

%% ===================== INDUCTANCE + BIOT-SAVART (pre-loop) ===============
[MArray, LTotal] = generate_inductance_matrix(numLines, numLinesAlongWire, ...
    wCond, thCond, x1Array, x2Array, y1Array, y2Array, z1Array, z2Array, ...
    mutualInductanceSpaceRatio, mutualInductanceSpaceRatioLimit, ...
    ignoreMutualInductancesTransverseElements);

MLPInv = create_MLPInv_matrix(MArray, numLines, numPoints, pointIndex1, pointIndex2);

% Initial current distribution: uniform along the wire, zero in shorts.
IArray = zeros(numLines, 1);
IArray(1:numLinesAlongWire) = I0;

[BiotSavartMatrixX, BiotSavartMatrixY, BiotSavartMatrixZ, ...
 BLocalXArray, BLocalYArray, BLocalZArray] = calc_biot_savart( ...
    numLines, x1Array, x2Array, y1Array, y2Array, z1Array, z2Array, ...
    IArray, mutualInductanceSpaceRatio);

%% ===================== INITIAL CONDITIONS ================================
temperatureArray = ones(numPoints - 1, numThermalSubdivisions) * initial_temperature;
% Initial hotspot, used to seed the quench:
temperatureArray(floor(numPoints/2 + 0.5), floor(numThermalSubdivisions/2 + 1.5)) = 22;

V_local_element = zeros(numLinesAlongWire, 1);
VGroundVector   = zeros(numLinesAlongWire, 1);

% Geometry-derived constants used in the heat-balance loop
arclengthCylinder  = 2*pi*RSol;
widthWall          = 0.005;     % m

A_helium  = (wCond - heightSlot) .* lenArrayAlongWireElements;
L_helium  = 0.5 .* ((wCond - heightSlot) + lenArrayAlongWireElements);
Nu        = 4.36;   % Nusselt number, steady laminar pipe flow

volumeArray = (lenArray(1:numLinesAlongWire) ./ numThermalSubdivisions ...
               .* thCond * wCond) * ones(1, numThermalSubdivisions);
Acylinder = 2 * 2*pi*RSol * (wCond - heightSlot) * numWindings;
Aslots    = 2 * numTransverse * thSlot * thWallatShort;
Atotal    = Acylinder + Aslots;
cooling_area_correction_factor = Atotal / Acylinder;

IExt = I0;   % External current (held constant in this run)

%% ===================== OUTPUT FILES (logfile + video) ====================
FileName = append(Description, ...
    string(datetime('now', 'TimeZone', 'local', 'Format', 'd-MMM-y')), '_', ...
    string(datetime('now', 'TimeZone', 'local', 'Format', 'HH.mm.ss')), '.txt');

fileID = -1;
if writedata == 1
    fileID = fopen(FileName, 'w');   % TODO: parameterize output directory
    write_log_header(fileID, I0, RSol, numWindings, N_shorts, resolution, ...
        thCond, HTStapeName, N_tapes, initial_temperature, PHeater, theater, ...
        dTMax, temperatureSubSteps, Helium_cooling, copper_heat, copper_conduction);
end

v = [];
if writevideo == 1
    v = VideoWriter(append(Description, '.mp4'), 'MPEG-4');
    v.FrameRate = 5;
    open(v);
end

%% ===================== TRANSIENT MAIN LOOP ===============================
disp('Transient calculation');

% Per-run pre-allocations.
% NOTE: many of these grow with iterationIndex; we pre-allocate to the cap
% (maxIteration) and the simulation may exit early via the time-limit break.
tArrayhistory                = zeros(1, maxIteration);
IArrayhistory                = zeros(numLines, maxIteration);
IcArrayhistory               = zeros(numPoints - 1, maxIteration);
VGroundVectorhistory         = zeros(numPoints - 1, maxIteration);
Qheliumcoolingarrayhistory   = zeros(numPoints - 1, numThermalSubdivisions, maxIteration);
temperatureArraystored       = zeros(numPoints - 1, numThermalSubdivisions, maxIteration);
temperatureRingTophistory    = zeros(1, maxIteration);
temperatureRingBottomhistory = zeros(1, maxIteration);
temperaturemax               = zeros(1, maxIteration);
centerBhistory               = zeros(1, maxIteration);
Efield                       = zeros(1, maxIteration);

s.MLPInv = MLPInv;

% Timing accumulators (use posixtime-equivalent: now*86400 in seconds)
overallFirstStartTime = now * 86400;

centerB = drawFieldMap(RSol, len, numPoints, numLines, numLinesAlongWire, ...
    numTransverse, numThermalSubdivisions, xPlot, yPlot, zPlot, ...
    x1Array, x2Array, y1Array, y2Array, z1Array, z2Array, ...
    BLocalXArray, BLocalYArray, BLocalZArray, ...
    temperatureArray, mutualInductanceSpaceRatio, IArray);

t                     = 0;                  % Initial simulated time (s)
temperatureRingTop    = temperatureRing;
temperatureRingBottom = temperatureRing;

% E-field probe arc length: distance between voltage taps at 1/10 and 9/10
% of the wire, projected onto the solenoid circumference.
EfieldArcLength = numWindings * arclengthCylinder * 0.8;
% NOTE: original hardcoded 4*2*pi*0.12 (= 4 turns * circumference at R=0.12).
% Replaced with EfieldArcLength = 0.8 * numWindings * 2*pi*RSol so the metric
% scales with the geometry; for the published 5-turn 0.12 m run this equals
% 4 * 2*pi*0.12 to within rounding.

for iterationIndex = 1:maxIteration
    overallStartTime              = now * 86400;
    materialsPropsCalculationStartTime = overallStartTime;

    externaldIdtVector             = zeros(numPoints - 1, 1);
    s.externalVoltagedIdtVector    = zeros(numPoints - 1, 1);

    %% --- Material property lookup ---------------------------------------
    % Index into the 0.1-K-spaced lookup tables.
    temperatureIndices         = round(10 * temperatureArray);
    thermalConductivityArray   = ktable(temperatureIndices);
    resistivityNormalArray     = rhotable(temperatureIndices);
    heatCapacityArray          = heat_capacity_al_alloy_5083(temperatureArray, volumeArray);

    BXYZarray   = [BLocalXArray, BLocalYArray, BLocalZArray];
    XYZarray    = [xArray, yArray, zArray];   %#ok<NASGU> kept for debugger
    BArray      = vecnorm(BXYZarray, 2, 2);

    % Tape c-axis is normal to the solenoid surface (z-component zeroed).
    CaxistapeArray   = [xArray, yArray, zeros(length(xArray), 1)];
    nC               = length(CaxistapeArray);
    theta_angleArray = acos(dot(BXYZarray(1:nC, :), CaxistapeArray, 2) ./ ...
                            (BArray(1:nC) .* vecnorm(CaxistapeArray, 2, 2)));

    % Worst-case (hottest) temperature in each line element drives Ic.
    maxTemparray = max(temperatureArray, [], 2);
    nT = length(temperatureArray);
    IcArray = N_tapes * tapewidth*1e3 * SCthickness*1e3 * ...
              parametrisation_fujikura(BArray(1:nT), maxTemparray, theta_angleArray(1:nT));
    IcArrayhistory(:, iterationIndex) = IcArray;

    materialsPropsCalculationFinishTime = now * 86400;

    %% --- Resistance of axial shorts (turn-to-turn) ----------------------
    temperatureBottomArray = temperatureArray(1:resolution:(numTransverse-1)*resolution, 1);
    temperatureTopArray    = temperatureArray( ...
        numPointsPerTurn+1 : resolution : (numTransverse-1)*resolution + numPointsPerTurn, 1);
    temperatureshort       = (temperatureBottomArray + temperatureTopArray) / 2;
    axialshortResistance   = heightSlot * rhotable(round(temperatureshort*10)) ...
                             / (thSlot * widthWall);
    axialshortResistance   = transpose(axialshortResistance);

    RArrayTrans          = axialshortResistance;
    RArrayTrans(end + 1) = RArrayTrans(end);   % padding to match numTransverse

    %% --- Longitudinal stabiliser resistance -----------------------------
    RNormalArray = resistivityNormalArray .* lenArrayAlongWire(1:numLinesAlongWire) ...
                   ./ (ACond .* numThermalSubdivisions);

    if copper_conduction == 1
        % First and last 1/5 of the conductor length are stabilised by the
        % copper rings (current leads). Override the local resistance.
        nFifth = round(length(RNormalArray) / 5);
        rhoCu_top    = rhoCu_nist(temperatureRingBottom*ones(nFifth, 1), ...
                                  BArray(1:nFifth), RingRRR, 1);
        rhoCu_bottom = rhoCu_nist(temperatureRingBottom*ones(nFifth, 1), ...
                                  BArray(1+4*nFifth:5*nFifth), RingRRR, 1);

        RNormalArray(1:nFifth, :) = ones(nFifth, 3) .* rhoCu_top ...
            .* lenArrayAlongWire(1:numLinesAlongWire/5) ...
            ./ (crossSectionRing .* numThermalSubdivisions);
        RNormalArray(1 + 4*nFifth : end, :) = ones(nFifth, 3) .* rhoCu_bottom ...
            .* lenArrayAlongWire(1:numLinesAlongWire/5) ...
            ./ (crossSectionRing .* numThermalSubdivisions);
    end

    %% --- Power dissipation per element (critical-state model) -----------
    IAbsArray    = abs(IArray(1:numLinesAlongWire)) * ones(1, numThermalSubdivisions);
    INormalArray = IAbsArray - IcArray;
    INormalArray = INormalArray .* (INormalArray > 0);
    VElementArray = INormalArray .* RNormalArray;
    PElementArray = VElementArray .* IAbsArray;

    % Joule heating in axial shorts split 50/50 between top and bottom turn.
    IArrayTrans              = IArray(numLinesAlongWire+1:end-1);
    PElementArrayTransverse  = IArrayTrans.^2 .* axialshortResistance;
    bottomIdx = 1 : resolution : (numTransverse-1)*resolution;
    topIdx    = numPointsPerTurn+1 : resolution : (numTransverse-1)*resolution + numPointsPerTurn;
    PElementArray(bottomIdx, 1) = PElementArray(bottomIdx, 1) + 0.5 * PElementArrayTransverse;
    PElementArray(topIdx,    1) = PElementArray(topIdx,    1) + 0.5 * PElementArrayTransverse;

    %% --- External heater pulses ----------------------------------------
    halfindex   = round(length(PElementArray) / 2);
    [PElementArray, Pheatersave] = apply_heater_pulses( ...
        PElementArray, halfindex, t, toffset, theater, PHeater);

    if t > simulation_time_limit
        break;
    end

    %% --- Thermal substep loop ------------------------------------------
    temperatureCalculationStartTime = now * 86400;
    dtInt = 0;

    for temperatureIterationIndex = 1:temperatureSubSteps
        % Right- and left-shifted neighbours for axial heat-flux stencil.
        temperatureArrayRight = [temperatureArray(:, 2:numThermalSubdivisions), ...
                                 temperatureArray(:, numThermalSubdivisions)];
        temperatureArrayRight(1:numPoints-2, numThermalSubdivisions) = ...
            temperatureArray(2:numPoints-1, 1);

        thermalConductivityArrayRight = [thermalConductivityArray(:, 2:numThermalSubdivisions), ...
                                         thermalConductivityArray(:, numThermalSubdivisions)];
        thermalConductivityArrayRight(1:numPoints-2, numThermalSubdivisions) = ...
            thermalConductivityArray(2:numPoints-1, 1);

        temperatureArrayLeft = [temperatureArray(:, 1), ...
                                temperatureArray(:, 1:numThermalSubdivisions-1)];
        temperatureArrayLeft(2:numPoints-1, 1) = ...
            temperatureArray(1:numPoints-2, numThermalSubdivisions);

        thermalConductivityArrayLeft = [thermalConductivityArray(:, 1), ...
                                        thermalConductivityArray(:, 1:numThermalSubdivisions-1)];
        thermalConductivityArrayLeft(2:numPoints-1, 1) = ...
            thermalConductivityArray(1:numPoints-2, numThermalSubdivisions);

        % Top/bottom temperatures and conductivities of axial shorts.
        temperatureBottomArray = temperatureArray(bottomIdx, 1);
        temperatureTopArray    = temperatureArray(topIdx,    1);
        thermalConductivityBottomArray = thermalConductivityArray(bottomIdx, 1);
        thermalConductivityTopArray    = thermalConductivityArray(topIdx,    1);

        % --- Axial conduction inside line element (right + left neighbour)
        QElementArray_temp = (temperatureArrayRight - temperatureArray) .* 0.5 ...
            .* (thermalConductivityArray + thermalConductivityArrayRight) ...
            .* ACond ./ lenArrayAlongWireElements;

        QElementArray_temp = QElementArray_temp + ...
            (temperatureArrayLeft - temperatureArray) .* 0.5 ...
            .* (thermalConductivityArray + thermalConductivityArrayLeft) ...
            .* ACond ./ lenArrayAlongWireElements;

        % --- Heat flow through axial shorts
        Qaxialarray = (temperatureTopArray - temperatureBottomArray) .* 0.5 ...
            .* (thermalConductivityBottomArray + thermalConductivityTopArray) ...
            .* thSlot .* thWallatShort ./ heightSlot;

        % --- Heat flow into the copper current leads
        % Approximation: local Al thermal conductivity is the limiter; copper
        % temperature is treated as uniform. Neglects Al/Cu contact resistance.
        QringBottomarray = (temperatureArray(1:numPointsPerTurn, :) - temperatureRingBottom) ...
            .* thermalConductivityArray(1:numPointsPerTurn, :) ...
            .* lenArrayAlongWireElements(1:numPointsPerTurn, :) * thCond ./ (0.5 * wCond);

        QringToparray = (temperatureArray(end+1-numPointsPerTurn:end, :) - temperatureRingTop) ...
            .* thermalConductivityArray(end+1-numPointsPerTurn:end, :) ...
            .* lenArrayAlongWireElements(end+1-numPointsPerTurn:end, :) * thCond ./ (0.5 * wCond);

        % --- Helium gas convective cooling
        Qheliumcoolingarray = zeros(size(temperatureArray));
        if Helium_cooling == 1
            k_heliumarray  = k_helium * ones(size(temperatureArray));
            h_helium_array = Nu * k_heliumarray ./ L_helium;
            R_thermal      = 0.5 * thCond ./ thermalConductivityArray + 1 ./ h_helium_array;
            Qheliumcoolingarray = 2 .* cooling_area_correction_factor .* A_helium ...
                                  .* (temperatureArray - temperature_helium) ./ R_thermal;
        end

        % --- Add heat flow from axial shorts (sign convention: + into bottom, - from top)
        QElementArray_temp(bottomIdx, 1) = QElementArray_temp(bottomIdx, 1) + Qaxialarray;
        QElementArray_temp(topIdx,    1) = QElementArray_temp(topIdx,    1) - Qaxialarray;

        % --- Add heat flow into copper rings
        if copper_heat == 1
            QElementArray_temp(1:numPointsPerTurn, :) = ...
                QElementArray_temp(1:numPointsPerTurn, :) - QringBottomarray;
            QElementArray_temp(end+1-numPointsPerTurn:end, :) = ...
                QElementArray_temp(end+1-numPointsPerTurn:end, :) - QringToparray;
        end

        QElementArray = QElementArray_temp - Qheliumcoolingarray;

        % NOTE: heat-capacity is already mass*cp (volume incorporated upstream),
        % so we divide by heatCapacityArray (not heatCapacityArray.*volumeArray).
        dTdtElementArray = (PElementArray + QElementArray) ./ heatCapacityArray;

        dTdtMax = max(max(abs(dTdtElementArray)));
        dtLocal = dTMax ./ dTdtMax / temperatureSubSteps;

        if dTdtMax / temperatureSubSteps == 0
            warning('main_5turnMM:zero_dTdt', ...
                'dT/dt is zero everywhere; thermal substep timestep is undefined.');
        end

        heatCapacityCoppertop    = CpCu_nist(temperatureRingTop,    massRing);
        heatCapacityCopperbottom = CpCu_nist(temperatureRingBottom, massRing);
        dTdtCoppertop    = sum(sum(QringToparray))    / heatCapacityCoppertop;
        dTdtCopperBottom = sum(sum(QringBottomarray)) / heatCapacityCopperbottom;

        dtInt                 = dtInt + dtLocal;
        temperatureArray      = temperatureArray + dtLocal * dTdtElementArray;
        temperatureRingBottom = temperatureRingBottom + dtLocal * dTdtCopperBottom;
        temperatureRingTop    = temperatureRingTop    + dtLocal * dTdtCoppertop;
    end

    temperatureArraystored(:, :, iterationIndex) = temperatureArray;
    dt = dtInt;
    temperatureRingTophistory(iterationIndex)    = temperatureRingTop;
    temperatureRingBottomhistory(iterationIndex) = temperatureRingBottom;
    temperaturemax(iterationIndex)               = max(temperatureArray(:));

    temperatureCalculationFinishTime = now * 86400;
    ODECalculationStartTime          = now * 86400;

    %% --- ODE for current redistribution --------------------------------
    s.RArrayTrans            = RArrayTrans;
    s.IcArray                = IcArray;
    s.bVector                = [];
    s.VGroundVector          = [];
    s.numLines               = numLines;
    s.numLinesAlongWire      = numLinesAlongWire;
    s.numTransverse          = numTransverse;
    s.numPoints              = numPoints;
    s.lenArrayAlongWire      = lenArrayAlongWire;
    s.signArray              = zeros(numLinesAlongWire, 1);
    s.INormalArray           = zeros(numLinesAlongWire, 1);
    s.RNormalArray           = RNormalArray;
    s.numThermalSubdivisions = numThermalSubdivisions;

    IArray0 = IArray;
    tspan   = [t, t + dt];
    myfun   = @(tt, II) myODE(tt, II, s);
    [tArchive, IArrayArchive] = ode15s(myfun, tspan, IArray0);
    t       = tArchive(end);
    IArray  = IArrayArchive(end, :)';

    tArrayhistory(iterationIndex)    = t;
    IArrayhistory(:, iterationIndex) = IArray;
    VGroundVector                    = myODEVoltages(t, IArray, s);

    Qheliumcoolingarrayhistory(:, :, iterationIndex) = Qheliumcoolingarray;
    VGroundVectorhistory(:, iterationIndex)          = VGroundVector;

    % E-field probe between 1/10 and 9/10 voltage taps.
    idxLow  = round(length(VGroundVector) * 1/10);
    idxHigh = round(length(VGroundVector) * 9/10);
    Efield(iterationIndex) = 1e6 * (VGroundVector(idxLow) - VGroundVector(idxHigh)) ...
                             / EfieldArcLength;

    ODECalculationFinishTime = now * 86400;

    %% --- Field-map update + visualisation ------------------------------
    plottingStartTime = now * 86400;

    BLocalXArray = BiotSavartMatrixX * IArray;
    BLocalYArray = BiotSavartMatrixY * IArray;
    BLocalZArray = BiotSavartMatrixZ * IArray;

    if mod(iterationIndex + 1, drawFigureAtIteration) == 0
        centerB = drawFieldMap(RSol, len, numPoints, numLines, numLinesAlongWire, ...
            numTransverse, numThermalSubdivisions, xPlot, yPlot, zPlot, ...
            x1Array, x2Array, y1Array, y2Array, z1Array, z2Array, ...
            BLocalXArray, BLocalYArray, BLocalZArray, ...
            temperatureArray, mutualInductanceSpaceRatio, IArray);
        title(sprintf('Iteration = %0.0f, t = %0.3g s, I_{ext} = %0.2f A', ...
            iterationIndex, t, IExt));
        pause(0.1);

        if writevideo == 1
            figure(1); hold on;
            frame = getframe(gcf);
            writeVideo(v, frame);
            hold off;
        end

        centerBhistory(iterationIndex) = centerB;
        figure(2);
        plot(tArrayhistory(centerBhistory > 0), centerBhistory(centerBhistory > 0), ...
            Color='b', Marker='.', MarkerSize=10, LineStyle='-');
        xlabel('t [s]'); ylabel('B_{center} [T]'); grid on;

        figure(3);
        plot(tArrayhistory, Efield, Color='b', Marker='.', MarkerSize=10, LineStyle='-');
        xlabel('t [s]'); ylabel('E_{field} [\muV/m]'); grid on;
    end

    plottingFinishTime = now * 86400;
    overallFinishTime  = now * 86400;

    %% --- Per-iteration timing report -----------------------------------
    deltaTotal  = max(overallFinishTime - overallStartTime, eps);
    fracMatCalc = 100 * (materialsPropsCalculationFinishTime - materialsPropsCalculationStartTime) / deltaTotal;
    fracTSteps  = 100 * (temperatureCalculationFinishTime  - temperatureCalculationStartTime)  / deltaTotal;
    fracODE     = 100 * (ODECalculationFinishTime - ODECalculationStartTime) / deltaTotal;
    fracPlot    = 100 * (plottingFinishTime       - plottingStartTime)       / deltaTotal;

    fprintf(['Index: %0.0f, time expired: %0.1f s, tsim %0.6g s, ', ...
             'last iteration %0.1f s ', ...
             '(%0.2f%% mat, %0.2f%% TSteps, %0.2f%% ODE, %0.2f%% plot)\n'], ...
        iterationIndex, overallFinishTime - overallFirstStartTime, t, ...
        deltaTotal, fracMatCalc, fracTSteps, fracODE, fracPlot);

    if writedata == 1
        fprintf(fileID, '%0.6g, %0.6f, %0.3f, %0.3f, %0.3f, %0.6f, %0.4f \n', ...
            [tArrayhistory(end), centerB, Pheatersave, ...
             temperatureRingTop, temperatureRingBottom, Efield(end), ...
             temperaturemax(end)]);
    end
end

%% ===================== CLEANUP / FINAL SAVE ==============================
if writedata == 1 && fileID > 0
    fclose(fileID);
end
if writevideo == 1 && ~isempty(v)
    close(v);
end

savename = append(Description, ...
    string(datetime('now', 'TimeZone', 'local', 'Format', 'd-MMM-y')), '_', ...
    string(datetime('now', 'TimeZone', 'local', 'Format', 'HH.mm.ss')), '.mat');
save(savename);   % TODO: parameterize output directory

toc;

%% ===================== LOCAL FUNCTIONS ===================================

function [PElementArray, Pheatersave] = apply_heater_pulses( ...
        PElementArray, halfindex, t, toffset, theater, PHeater)
    % APPLY_HEATER_PULSES Add forced heater power to the central element
    % during each scheduled pulse window. Pulse k is active during the
    % interval [tStart_k, tStart_k + theater(k)], where tStart_k is the
    % cumulative offset built from toffset/theater (see definition below).
    Pheatersave = 0;

    pulseStarts    = zeros(1, numel(theater));
    pulseStarts(1) = toffset(1);
    for k = 2:numel(theater)
        pulseStarts(k) = pulseStarts(k-1) + theater(k-1) + toffset(k);
    end

    for k = 1:numel(theater)
        if pulseStarts(k) < t && t < pulseStarts(k) + theater(k)
            % Heater power is split across all numThermalSubdivisions of the
            % central line element. The "/3" assumes numThermalSubdivisions==3.
            PElementArray(halfindex, :) = PElementArray(halfindex, :) + PHeater / 3;
            Pheatersave = PHeater;
            return;
        end
    end
end

function write_log_header(fileID, I0, RSol, numWindings, N_shorts, resolution, ...
        thCond, HTStapeName, N_tapes, initial_temperature, PHeater, theater, ...
        dTMax, temperatureSubSteps, Helium_cooling, copper_heat, copper_conduction)
    % WRITE_LOG_HEADER Emit a fixed-format header to the run log.
    fprintf(fileID, 'Input Parameters RePISoSi code:\n\n');
    fprintf(fileID, 'Initial Current: %.0f [A]\n',                       I0);
    fprintf(fileID, 'Radius solenoid: %.3f [m]\n',                       RSol);
    fprintf(fileID, 'Number of turns: %.1g [-]\n',                       numWindings);
    fprintf(fileID, 'Number of shorts per turn: %.2g [-]\n',             N_shorts);
    fprintf(fileID, 'Number of line-elements between axial shorts: %.1g [-]\n', resolution);
    fprintf(fileID, 'Thickness of stabiliser wall: %.2g [mm]\n',         1e3*thCond);
    fprintf(fileID, 'Tape type: %s\n',                                   HTStapeName);
    fprintf(fileID, 'Number of tapes [-]: %.1g\n',                       N_tapes);
    fprintf(fileID, 'Initial temperature [K]: %.1f\n',                   initial_temperature);
    fprintf(fileID, 'Power of Heater [W]: %.1f\n',                       PHeater);
    fprintf(fileID, 'Heating time [s]: %.1f\n',                          theater);
    fprintf(fileID, 'deltaTmax (timestep determination) [K]: %.1f\n',    dTMax);
    fprintf(fileID, 'temperatureSubSteps [-]: %.1f\n',                   temperatureSubSteps);
    fprintf(fileID, 'Status Helium convection cooling (1=on,0=off): %.1g\n',     Helium_cooling);
    fprintf(fileID, 'Status Copper heat conduction cooling (1=on,0=off): %.1g\n', copper_heat);
    fprintf(fileID, 'Status Copper conductivitiy (1=on,0=off): %.1g\n\n\n',       copper_conduction);
    fprintf(fileID, 'time [s], centerB [T], Pheater [W], TringTop [K], TringBottom [K], Efield [uV/m], Tmax [K]\n');
end
