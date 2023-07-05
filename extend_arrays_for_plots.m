function [xPlot, yPlot, zPlot] = extend_arrays_for_plots(numPoints, ... 
    numThermalSubdivisions, xArray, yArray, zArray, x1Array, x2Array, y1Array, y2Array, z1Array, z2Array, wCond, thCond, RSolEff) 

    xArrayExtended = zeros((numPoints - 1)*numThermalSubdivisions + 1, 1);
    yArrayExtended = zeros((numPoints - 1)*numThermalSubdivisions + 1, 1);
    zArrayExtended = zeros((numPoints - 1)*numThermalSubdivisions + 1, 1);

    for index1 = 1:numPoints-1
        xArrayLocal = linspace(x1Array(index1), x2Array(index1), numThermalSubdivisions + 1);
        yArrayLocal = linspace(y1Array(index1), y2Array(index1), numThermalSubdivisions + 1);
        zArrayLocal = linspace(z1Array(index1), z2Array(index1), numThermalSubdivisions + 1);

        xArrayExtended((index1 - 1)*numThermalSubdivisions + (1:numThermalSubdivisions), 1) = xArrayLocal(1:numThermalSubdivisions);
        yArrayExtended((index1 - 1)*numThermalSubdivisions + (1:numThermalSubdivisions), 1) = yArrayLocal(1:numThermalSubdivisions);
        zArrayExtended((index1 - 1)*numThermalSubdivisions + (1:numThermalSubdivisions), 1) = zArrayLocal(1:numThermalSubdivisions);
    end

    xArrayExtended((numPoints - 1)*numThermalSubdivisions + 1) = xArray(numPoints);
    yArrayExtended((numPoints - 1)*numThermalSubdivisions + 1) = yArray(numPoints);
    zArrayExtended((numPoints - 1)*numThermalSubdivisions + 1) = zArray(numPoints);

    fmin = (RSolEff - 0.5*thCond)/RSolEff;
    fplus = (RSolEff + 0.5*thCond)/RSolEff;
    xPlot = [xArrayExtended*fmin, xArrayExtended*fmin, xArrayExtended*fplus, xArrayExtended*fplus, xArrayExtended*fmin];
    yPlot = [yArrayExtended*fmin, yArrayExtended*fmin, yArrayExtended*fplus, yArrayExtended*fplus, yArrayExtended*fmin];
    zPlot = [zArrayExtended - 0.5*wCond, zArrayExtended + 0.5*wCond, zArrayExtended + 0.5*wCond, zArrayExtended - 0.5*wCond, zArrayExtended - 0.5*wCond];

end 