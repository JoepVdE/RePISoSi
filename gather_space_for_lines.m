% -------------------------------------------------------------------------
%  gather_space_for_lines.m  --  Build the (x1, x2, y1, y2, z1, z2) end-
%                                point arrays for every line element in the
%                                solenoid mesh, including off-by-one-turn
%                                axial shorts.
%  Part of RePISoSi - https://github.com/JoepVdE/RePISoSi  -  License: MIT
%  Author : J.L. Van den Eijnden, 2026
%           Original implementation by M. Mentink (Oct. 2022),
%           extended by J.L. Van den Eijnden and A. Vaskuri (Oct. 2023).
%
%  Inputs:
%      numLinesAlongWire : number of straight line segments tracking the
%                          conductor helix.
%      numTransverse     : number of axial shorts (turn-to-turn).
%      numLines          : numLinesAlongWire + numTransverse.
%      numPointsPerTurn  : number of mesh points per helix turn.
%      xArray, yArray, zArray : Cartesian coordinates of the helix points
%                               (column vectors of length numPoints).
%      resolution        : number of unconnected line segments between two
%                          consecutive axial shorts (must be even).
%
%  Outputs:
%      x1Array, ..., z2Array : end-point arrays (numLines x 1 each).
%      pointIndex1, pointIndex2 : 1-based indices into the helix point list
%                                 identifying the two endpoints of every line.
%
%  Each axial short connects point (resolution*(k-1)+1) to point
%  (resolution*(k-1)+1 + numPointsPerTurn). The "+ numPointsPerTurn" gives
%  the off-by-one-turn pattern that maximises axial resistance.
% -------------------------------------------------------------------------
function [x1Array, x2Array, y1Array, y2Array, z1Array, z2Array, ...
          pointIndex1, pointIndex2] = gather_space_for_lines( ...
            numLinesAlongWire, numTransverse, numLines, numPointsPerTurn, ...
            xArray, yArray, zArray, resolution)

    % Pre-allocate
    x1Array     = zeros(numLines, 1);
    x2Array     = zeros(numLines, 1);
    y1Array     = zeros(numLines, 1);
    y2Array     = zeros(numLines, 1);
    z1Array     = zeros(numLines, 1);
    z2Array     = zeros(numLines, 1);
    pointIndex1 = zeros(numLines, 1);
    pointIndex2 = zeros(numLines, 1);

    %% --- Lines tracking the conductor (helix) --------------------------
    for index = 1:numLinesAlongWire
        x1Array(index) = xArray(index);
        x2Array(index) = xArray(index + 1);
        y1Array(index) = yArray(index);
        y2Array(index) = yArray(index + 1);
        z1Array(index) = zArray(index);
        z2Array(index) = zArray(index + 1);

        pointIndex1(index) = index;
        pointIndex2(index) = index + 1;
    end

    %% --- Axial shorts (off-by-one-turn) --------------------------------
    spacing = resolution;   % unconnected line segments between two shorts
    for index = 1:numTransverse
        indexRef  = index + numLinesAlongWire;
        startIdx  = spacing*(index - 1) + 1;
        endIdx    = startIdx + numPointsPerTurn;   % off-by-one-turn

        x1Array(indexRef) = xArray(startIdx);
        x2Array(indexRef) = xArray(endIdx);
        y1Array(indexRef) = yArray(startIdx);
        y2Array(indexRef) = yArray(endIdx);
        z1Array(indexRef) = zArray(startIdx);
        z2Array(indexRef) = zArray(endIdx);

        pointIndex1(indexRef) = startIdx;
        pointIndex2(indexRef) = startIdx + numPointsPerTurn;   % off-by-one-turn
    end
end
