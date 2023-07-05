function [BPoints] = calcPoints(numFieldPoints, numLines, xPointArray, yPointArray, zPointArray, ... 
    xMidArray, yMidArray, zMidArray, xDiffArray, yDiffArray, zDiffArray, lenArray, IArray, ... 
    mutualInductanceSpaceRatio) 

    BiotSavartPointMatrix = zeros(numFieldPoints, numLines, 3);
    BPoints = zeros(numFieldPoints, 3);
    BPoints(:, 1) = xPointArray;
    BPoints(:, 2) = yPointArray;
    BPoints(:, 3) = zPointArray;

    % Here index denotes the element on which the field is calculated, whereas
    % index2 indicates the element producing the field
    % Note that self-field is ignored
    for index = 1:numFieldPoints
        r = sqrt((xPointArray(index)-xMidArray).^2+(yPointArray(index)-yMidArray).^2+(zPointArray(index)-zMidArray).^2);
        rInv3 = 1./r.^3;
        ratio = r./lenArray;

        for index2 = 1:numFieldPoints
            if(r(index2) == 0)
                rInv3(index2) = 0;
                ratio(index2) = mutualInductanceSpaceRatio;
            end        
            if(ratio(index2) < mutualInductanceSpaceRatio)
                numSubdivisions = floor(mutualInductanceSpaceRatio./ratio(index2)) + 1;
                x1 = xPointArray(index);
                y1 = yPointArray(index);
                z1 = zPointArray(index);

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

        BMag = 1e-7.*rInv3;
        v1 = [xDiffArray yDiffArray zDiffArray];
        v2 = [xPointArray(index) - xMidArray, yPointArray(index) - yMidArray, zPointArray(index) - zMidArray];
        BiotSavartPointMatrix(index, :, :) = (BMag*ones(1, 3)) .* cross(v1, v2);    
    end

    BiotSavartPointMatrixX = BiotSavartPointMatrix(:, :, 1);
    BiotSavartPointMatrixY = BiotSavartPointMatrix(:, :, 2);
    BiotSavartPointMatrixZ = BiotSavartPointMatrix(:, :, 3);

    BLocalXPointArray = BiotSavartPointMatrixX*IArray;
    BLocalYPointArray = BiotSavartPointMatrixY*IArray;
    BLocalZPointArray = BiotSavartPointMatrixZ*IArray;

    BPoints = [BPoints, BLocalXPointArray, BLocalYPointArray, BLocalZPointArray];
end 
