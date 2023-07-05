function [BiotSavartMatrixX, BiotSavartMatrixY, BiotSavartMatrixZ, BLocalXArray, BLocalYArray, BLocalZArray] ... 
    = calc_biot_savart(numLines, x1Array, x2Array, y1Array, y2Array, z1Array, z2Array, IArray, mutualInductanceSpaceRatio)
    % ***** Biot-Savart calculations ***** 
    % The special case of a uniform constant current I (current is out from the integral)
    disp("Biot-Savart Geometry calculation");
    f = waitbar(0, 'Biot-Savart Geometry calculation');
    
    updateCounter = floor(numLines/100);

    BiotSavartMatrix = zeros(numLines, numLines, 3); % Initialize matrix with zeros 

    xMidArray = 0.5*(x1Array + x2Array);
    yMidArray = 0.5*(y1Array + y2Array);
    zMidArray = 0.5*(z1Array + z2Array);

    % Distance between the turns             
    xDiffArray = x2Array - x1Array;
    yDiffArray = y2Array - y1Array;
    zDiffArray = z2Array - z1Array; 

    lenArray = sqrt(xDiffArray.^2 + yDiffArray.^2 + zDiffArray.^2);
    
    % Index denotes the element on which the field is calculated 
    % Index2 denotes the element producing the field
    % Self-field of the conductor is ignored
    for index = 1:numLines
        r = sqrt((xMidArray(index) - xMidArray).^2 + (yMidArray(index) - yMidArray).^2 + (zMidArray(index) - zMidArray).^2); % Distance from conductor r
        rInv3 = 1./(r.^3); % Inverse of the 3rd power of the distance from conductor 
        ratio = r./lenArray;

        for index2 = 1:numLines 
            if(r(index2) == 0)
                rInv3(index2) = 0;
                ratio(index2) = mutualInductanceSpaceRatio;
            end        
            if(ratio(index2) < mutualInductanceSpaceRatio)
                numSubdivisions = floor(mutualInductanceSpaceRatio./ratio(index2)) + 1;
                x1 = xMidArray(index);
                y1 = yMidArray(index);
                z1 = zMidArray(index);

                x20 = xMidArray(index2);
                y20 = yMidArray(index2);
                z20 = zMidArray(index2);            
                dx2 = xDiffArray(index2)/numSubdivisions;
                dy2 = yDiffArray(index2)/numSubdivisions;
                dz2 = zDiffArray(index2)/numSubdivisions;

                rInv3Local = 0;

                for index3 = 1:numSubdivisions
                    x2 = (index3 - 0.5)*dx2 + x20;
                    y2 = (index3 - 0.5)*dy2 + y20;
                    z2 = (index3 - 0.5)*dz2 + z20;

                    rLocal = sqrt((x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2);
                    rInv3Local = rInv3Local + 1/(rLocal.^3*numSubdivisions);
                end
                rInv3(index2) = rInv3Local;
            end
        end

        BMag = 1e-7.*rInv3; % mu0/4*pi = 1e-7 H/m
        v1 = [xDiffArray yDiffArray zDiffArray]; % Differential length dl per axis 
        v2 = [(xMidArray(index) - xMidArray), (yMidArray(index) - yMidArray), (zMidArray(index) - zMidArray)]; % Distance from conductor r per axis
        BiotSavartMatrix(index, :, :) = (BMag*ones(1, 3)) .* cross(v1, v2); % Cross product 

        if(mod(index,updateCounter) == 0)
            waitbar(index/numLines, f, sprintf("Biot-Savart Geometry calculation, %0.0f/100", index/numLines*100));
            pause(0.1);
        end
    end

    BiotSavartMatrixX = BiotSavartMatrix(:, :, 1); 
    BiotSavartMatrixY = BiotSavartMatrix(:, :, 2); 
    BiotSavartMatrixZ = BiotSavartMatrix(:, :, 3); 
    close(f); 

    BLocalXArray = BiotSavartMatrixX*IArray; % Magnetic flux density vector along x-axis (T)
    BLocalYArray = BiotSavartMatrixY*IArray; % Magnetic flux density vector along y-axis (T)
    BLocalZArray = BiotSavartMatrixZ*IArray; % Magnetic flux density vector along z-axis (T)
    % ***** 

    
    
end 