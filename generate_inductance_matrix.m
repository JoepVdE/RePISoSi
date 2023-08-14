function [MArray, LTotal] = generate_inductance_matrix(numLines, numLinesAlongWire, wCond, thCond, ... 
    x1Array, x2Array, y1Array, y2Array, z1Array, z2Array, ... 
    mutualInductanceSpaceRatio, mutualInductanceSpaceRatioLimit, ignoreMutualInductancesTransverseElements) 
    mu0 = pi*4E-7; %[H/m]


    xMidArray = 0.5*(x1Array + x2Array);
    yMidArray = 0.5*(y1Array + y2Array);
    zMidArray = 0.5*(z1Array + z2Array);

    % Distance between the turns             
    xDiffArray = x2Array - x1Array;
    yDiffArray = y2Array - y1Array;
    zDiffArray = z2Array - z1Array; 
    
    lenArray = sqrt(xDiffArray.^2 + yDiffArray.^2 + zDiffArray.^2);
    
    MArray = zeros(numLines, numLines);

    % ***** Generate inductance matrix ***** 
    f = waitbar(0, 'Generating inductance matrix');
    disp("Generating inductance matrix");

    updateCounter = floor(numLines/100);

    for index = 1:numLines 
        for index2 = 1:numLines
            if(index == index2) 
                % Self inductance
                MArray(index, index2) = calcSelfInductanceRectangularBar(wCond, thCond, lenArray(index));
            else 
                % Mutual inductance
                rInitial = sqrt((xMidArray(index) - xMidArray(index2)).^2 + (yMidArray(index) - yMidArray(index2)).^2 + (zMidArray(index) - zMidArray(index2)).^2);
                rInv = 1./rInitial;

                r1 = sqrt((x1Array(index) - x2Array(index2)).^2 + (y1Array(index) - y2Array(index2)).^2 + (z1Array(index) - z2Array(index2)).^2);
                r2 = sqrt((x2Array(index) - x1Array(index2)).^2 + (y2Array(index) - y1Array(index2)).^2 + (z2Array(index) - z1Array(index2)).^2);

                rMin = min([rInitial r1 r2]);

                lenAvg = 0.5*(lenArray(index) + lenArray(index2));

                if(rMin/lenAvg < mutualInductanceSpaceRatio)

                    if(rMin == 0)
                        numSubdivisions = mutualInductanceSpaceRatioLimit;
                    else
                        numSubdivisions = floor(mutualInductanceSpaceRatio/(rMin/lenAvg)) + 1;
                    end
                    rInv = 0; 

                    x10 = x1Array(index);
                    y10 = y1Array(index);
                    z10 = z1Array(index);
                    x20 = x1Array(index2);
                    y20 = y1Array(index2);
                    z20 = z1Array(index2);                

                    dx1 = xDiffArray(index)/numSubdivisions;
                    dy1 = yDiffArray(index)/numSubdivisions;
                    dz1 = zDiffArray(index)/numSubdivisions;
                    dx2 = xDiffArray(index2)/numSubdivisions;
                    dy2 = yDiffArray(index2)/numSubdivisions;
                    dz2 = zDiffArray(index2)/numSubdivisions;

                    for index3 = 1:numSubdivisions
                        for index4 = 1:numSubdivisions
                            x1 = (index3 - 0.5)*dx1 + x10;
                            y1 = (index3 - 0.5)*dy1 + y10;
                            z1 = (index3 - 0.5)*dz1 + z10;
                            x2 = (index4 - 0.5)*dx2 + x20;
                            y2 = (index4 - 0.5)*dy2 + y20;
                            z2 = (index4 - 0.5)*dz2 + z20;

                            r = sqrt((x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2);
                            rInv = rInv + 1/(r*numSubdivisions^2);
                        end
                    end
                end

                v1 = [xDiffArray(index), yDiffArray(index), zDiffArray(index)];
                v2 = [xDiffArray(index2), yDiffArray(index2), zDiffArray(index2)]; 

                MArray(index, index2) = mu0*dot(v1, v2)*rInv/(4*pi); 
            end
        end

        % ***** Plot progress bar of inductance matrix calculation ***** 
        if(mod(index, updateCounter) == 0)
            waitbar(index/numLines, f, sprintf("Generating inductance matrix, %0.0f/100", index/numLines*100));
            pause(0.1);
        end
        % *****
    end


    if(ignoreMutualInductancesTransverseElements == 1)
        % Mutual inductance at the matrix diagonal
        MArray2 = zeros(numLines, numLines);
        for index1 = 1:numLines
            MArray2(index1, index1) = MArray(index1, index1);
        end

        for index1 = 1:numLinesAlongWire
            for index2 = 1:numLinesAlongWire
                MArray2(index1, index2) = MArray(index1, index2);
            end
        end
        MArray = MArray2;
    end

    % Mutual inductance checks
    kMax = 0;
    kMaxIndex1 = 0;
    kMaxIndex2 = 0;


    for index = 1:numLines
        for index2 = 1:numLines
            if(index ~= index2)
                %kLocal(index, index2) = MArray(index, index2)/sqrt(MArray(index, index)*MArray(index2, index2));
                kLocal = MArray(index, index2)/sqrt(MArray(index, index)*MArray(index2, index2));
                if(kLocal > kMax)
                    kMax = kLocal;
                    kMaxIndex1 = index;
                    kMaxIndex2 = index2;
                end
            end
        end
    end

    close(f);
    
    LTotal = sum(sum(MArray(1:numLinesAlongWire, 1:numLinesAlongWire)));
    fprintf("Inductance, excluding transverse elements: %0.5g\n", LTotal);
    fprintf("Inductance, including transverse elements: %0.5g\n", sum(sum(MArray)));
    fprintf("Highest coupling constant: %0.3g\n", kMax);
    if(kMax > 1)
        disp(["Incorrect coupling matrix, kMaxIndex1, kMaxIndex2",kMaxIndex1,kMaxIndex1]);
        return
    end 
    
end 
