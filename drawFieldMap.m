function [] = drawFieldMap(RSol, len, numPoints, numLines, numLinesAlongWire, numThermalSubdivisions, xPlot, yPlot, zPlot, ... 
    x1Array, x2Array, y1Array, y2Array, z1Array, z2Array, BLocalXArray, BLocalYArray, BLocalZArray, ... 
    temperatureArray, mutualInductanceSpaceRatio, IArray) 
    % drawFieldMap is a script relying on the parameters specified in run
    % script 

    % Distance between the turns             
    xDiffArray = x2Array - x1Array;
    yDiffArray = y2Array - y1Array;
    zDiffArray = z2Array - z1Array; 

    lenArray = sqrt(xDiffArray.^2 + yDiffArray.^2 + zDiffArray.^2);
%     lenArrayAlongWire = lenArray(1:numLinesAlongWire);
%     lenArrayAlongWireElements = lenArrayAlongWire*ones(1, numThermalSubdivisions)./numThermalSubdivisions;

    % Middle coordinates between the turns 
    xMidArray = 0.5*(x1Array + x2Array);
    yMidArray = 0.5*(y1Array + y2Array);
    zMidArray = 0.5*(z1Array + z2Array); 
    
    
    
    set(gcf, 'color', 'w'); 
    hold off;

    % Plot selections (Yes = 1, No = 0) 
    plotLongitudinalLines = 1;     
    plotLongitudinalLinesBMag = 0;  
    plotTransverseLines = 0;
    plotFieldOnConductor = 1;
    plotCurrent = 1;
    plotFieldMap = 1;

    % Color transparency
    conductorTransparency = 0.5;

    numPointsFieldMap = 21;

    % Scaling of the axes 
    maxSize = 1.5*max(RSol, len);

    tmArray = zeros((numPoints - 1)*numThermalSubdivisions, 1);	

    for index = 1:(numPoints - 1)
        tmArray((index - 1)*numThermalSubdivisions + (1:numThermalSubdivisions), 1) = temperatureArray(index, :);
    end    
    tmArrayE = [tmArray; tmArray((numPoints - 1)*numThermalSubdivisions)];
    colorPlot = [tmArrayE, tmArrayE, tmArrayE, tmArrayE, tmArrayE]; 

    % First the lines 
    if(plotLongitudinalLines == 1)
        hold on;
        %for index=1:numLinesAlongWire
        %    plot3([x1Array(index) x2Array(index)],[y1Array(index) y2Array(index)],[z1Array(index) z2Array(index)],'b');
        %    hold on;
        %end
        %plot3([x1Array x2Array(numLinesAlongWire)],[y1Array y2Array(numLinesAlongWire)],[z1Array z2Array(numLinesAlongWire)],'b');
        %plot3([x1Array(1:numLinesAlongWire); x2Array(numLinesAlongWire)],[y1Array(1:numLinesAlongWire); y2Array(numLinesAlongWire)],[z1Array(1:numLinesAlongWire); z2Array(numLinesAlongWire)],'k');
        colormap('jet')
        daspect([1, 1, 1]);
        view(80, 40);


        if(plotCurrent == 1)
            sPlot = surf(xPlot, yPlot, zPlot, colorPlot, 'FaceAlpha', conductorTransparency);
        else
            sPlot = surf(xPlot, yPlot, zPlot, colorPlot);
        end
        %sPlot.EdgeColor = 'none';
        colormap('jet')
        handle = colorbar;
    end

    % First the lines
    if(plotLongitudinalLinesBMag == 1)
        BNormArray = sqrt(BLocalXArray.^2 + BLocalYArray.^2 + BLocalZArray.^2);
        BNormArray = BNormArray(1:numLinesAlongWire);
        %for index = 1:numLinesAlongWire
        %    plot3([x1Array(index) x2Array(index)],[y1Array(index) y2Array(index)],[z1Array(index) z2Array(index)],'b');
        %    hold on;
        %end
        %plot3([x1Array x2Array(numLinesAlongWire)],[y1Array y2Array(numLinesAlongWire)],[z1Array z2Array(numLinesAlongWire)],'b');
        p = plot3([x1Array(1:numLinesAlongWire); x2Array(numLinesAlongWire)], ...
                  [y1Array(1:numLinesAlongWire); y2Array(numLinesAlongWire)], ... 
                  [z1Array(1:numLinesAlongWire); z2Array(numLinesAlongWire)], 'k');
        hold on;
        %scatter3([xMidArray(1:numLinesAlongWire)],[yMidArray(1:numLinesAlongWire)],[zMidArray(1:numLinesAlongWire)],100,BNormArray,'filled','MarkerEdgeColor','b');
        p.LineWidth = 1;
        %view(40, 35)
        %handle = colorbar; 
    end

    if(plotTransverseLines > 0)
        for index = (numLinesAlongWire + 1):(numLinesAlongWire + numTransverse)
            plot3([x1Array(index); x2Array(index)], [y1Array(index); y2Array(index)], [z1Array(index); z2Array(index)], 'b');
            %plot3([x1Array(index) x2Array(index)],[y1Array(index) y2Array(index)],[z1Array(index) z2Array(index)],'-o','Color','b','MarkerSize',4)
            hold on;
        end
    end



    % Then the lines showing the field direction and magnitude

    BNormArray = sqrt(BLocalXArray.^2 + BLocalYArray.^2 + BLocalZArray.^2);

    %lenLineAvg = sum(lenArray)/numLines;
    %BAvg = sum(BNormArray)/numLines;

    BNormFactor = min(lenArray)./max(BNormArray);

    if(plotFieldOnConductor > 0)
        for index = 1:numLinesAlongWire
            x = xMidArray(index);
            y = yMidArray(index);
            z = zMidArray(index);

            x1 = x - BLocalXArray(index)*BNormFactor;
            x2 = x + BLocalXArray(index)*BNormFactor;
            y1 = y - BLocalYArray(index)*BNormFactor;
            y2 = y + BLocalYArray(index)*BNormFactor;
            z1 = z - BLocalZArray(index)*BNormFactor;
            z2 = z + BLocalZArray(index)*BNormFactor;

            plot3([x1, x2], [y1, y2], [z1, z2], 'r');
        end
    end

    if(plotFieldMap > 0)
        xmin = 0; 
        xmax = 0; 
        numx = 1;
        ymin = -RSol*3; 
        ymax = RSol*3; 
        numy = numPointsFieldMap;    
        zmin = -len*2; 
        zmax = len*2; 
        numz = numPointsFieldMap;

        disp('Plotting field map')

        [x, y, z] = meshgrid(linspace(xmin, xmax, numx), linspace(ymin, ymax, numy), linspace(zmin, zmax, numz));

        u = zeros(numy, numx, numz);
        v = zeros(numy, numx, numz);
        w = zeros(numy, numx, numz);

        numFieldPoints = numz;

        for indexx = 1:numx
            for indexy = 1:numy
                xPointArray = squeeze(x(indexy, indexx, :));
                yPointArray = squeeze(y(indexy, indexx, :));
                zPointArray = squeeze(z(indexy, indexx, :));

                BPoints = calcPoints(numFieldPoints, numLines, xPointArray, yPointArray, zPointArray, ... 
                    xMidArray, yMidArray, zMidArray, xDiffArray, yDiffArray, zDiffArray, lenArray, IArray, ... 
                    mutualInductanceSpaceRatio); 



                u(indexy, indexx, :) = BPoints(:, 4);
                v(indexy, indexx, :) = BPoints(:, 5);
                w(indexy, indexx, :) = BPoints(:, 6);
            end
        end

        norm = sqrt(u.^2 + v.^2 + w.^2); 
        colormap(jet); 
        quiverC3D(x, y, z, u, v, w, 1); 
        xlabel('x'); 
        ylabel('y'); 
        zlabel('z'); 
    end 

    currentScalingFactor = 0.5; 
    if(plotCurrent > 0) 
        %IMax = max(abs(IArray)); 
        IMaxFactor = min(lenArray./abs(IArray)); 

        unitVector = [x2Array - x1Array, y2Array - y1Array, z2Array - z1Array] ...
            ./sqrt((x2Array - x1Array).^2 + (y2Array - y1Array).^2 + (z2Array - z1Array).^2);
        %unitX = unitVector(:, 1).*IArray./IMax; 
        %unitY = unitVector(:, 2).*IArray./IMax; 
        %unitZ = unitVector(:, 3).*IArray./IMax; 

        %unitX = unitVector(:, 1).*IArray; 
        %unitY = unitVector(:, 2).*IArray; 
        %unitZ = unitVector(:, 3).*IArray; 

        unitX = unitVector(:, 1).*IArray.*IMaxFactor; 
        unitY = unitVector(:, 2).*IArray.*IMaxFactor; 
        unitZ = unitVector(:, 3).*IArray.*IMaxFactor; 

        %quiverC3D(x1Array,y1Array,z1Array,unitX,unitY,unitZ,'LineWidth',1,'MaxHeadSize',10);
        hold off;
        %quiverC3D(x1Array,y1Array,z1Array,unitX,unitY,unitZ,'LineWidth',1,'MaxHeadSize',10,'scale',0);
        quiver3(x1Array, y1Array, z1Array, unitX, unitY, unitZ, 'g');
        hold on;

        sPlot = surf(xPlot, yPlot, zPlot, colorPlot, 'FaceAlpha', conductorTransparency);
        %sPlot.EdgeColor = 'none';
        colormap('jet')
        handle = colorbar; 
        daspect([1, 1, 1]);
        view(80, 40);
    end

    axis([-maxSize, maxSize, -maxSize, maxSize, -maxSize, maxSize]);
    xlabel('{\it x} / m');
    ylabel('{\it y} / m');
    zlabel('{\it z} / m');
    ylabel(handle, 'Temperature / K'); 

end 
