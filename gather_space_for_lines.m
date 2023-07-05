function [x1Array, x2Array, y1Array, y2Array, z1Array, z2Array, pointIndex1, pointIndex2] = gather_space_for_lines(numLinesAlongWire, numTransverse, numLines, numPointsPerTurn, xArray, yArray, zArray, resolution) 



% Initialize coordinates around the conductor 
x1Array = zeros(numLines, 1);
x2Array = zeros(numLines, 1);
y1Array = zeros(numLines, 1);
y2Array = zeros(numLines, 1);
z1Array = zeros(numLines, 1);
z2Array = zeros(numLines, 1);

% xTransArray = zeros(numLines, 1);
% yTransArray = zeros(numLines, 1);
% zTransArray = zeros(numLines, 1);

pointIndex1 = zeros(numLines, 1);
pointIndex2 = zeros(numLines, 1);

% ***** Draw lines ***** 
for index = 1:numLinesAlongWire
    x1Array(index) = xArray(index);
    x2Array(index) = xArray(index + 1);
    
    y1Array(index) = yArray(index);
    y2Array(index) = yArray(index + 1);
    
    z1Array(index) = zArray(index);
    z2Array(index) = zArray(index + 1);
    
    pointIndex1(index) = index;
    pointIndex2(index) = index + 1;
    
    diffX = x2Array(index) - x1Array(index);
    diffY = y2Array(index) - y1Array(index);
    
    %lenTrans = sqrt(diffX.^2 + diffY.^2);
    
%     xTransArray(index) = diffY/lenTrans;
%     yTransArray(index) = -diffX/lenTrans;
%     zTransArray(index) = 0; % This assumes that the tapes are always aligned perpendicular to z 
    
end

% for index = 1:numTransverse
%     indexRef = index + numLinesAlongWire; 
%     x1Array(indexRef) = xArray(index); 
%     x2Array(indexRef) = xArray(index + numPointsPerTurn);   % Off-by-one-turn
%     y1Array(indexRef) = yArray(index);
%     y2Array(indexRef) = yArray(index + numPointsPerTurn);   % Off-by-one-turn
%     z1Array(indexRef) = zArray(index);
%     z2Array(indexRef) = zArray(index + numPointsPerTurn);   % Off-by-one-turn
%     
%     pointIndex1(indexRef) = index;
%     pointIndex2(indexRef) = index + numPointsPerTurn;       % Off-by-one-turn
% end




% %%%working for main
% for index = 1:numTransverse
%     spacing = (Ratio_odd_integer+1);
%     %spacing =2;
%     indexRef = index + numLinesAlongWire;
%     x1Array(indexRef) = xArray(spacing*(index-1)+1); 
%     x2Array(indexRef) = xArray(spacing*(index-1)+1 + numPointsPerTurn);   % Off-by-one-turn (and 1 element)
%     y1Array(indexRef) = yArray(spacing*(index-1)+1);
%     y2Array(indexRef) = yArray(spacing*(index-1)+1 + numPointsPerTurn);   % Off-by-one-turn
%     z1Array(indexRef) = zArray(spacing*(index-1)+1);
%     z2Array(indexRef) = zArray(spacing*(index-1)+1 + numPointsPerTurn);   % Off-by-one-turn
% 
%     pointIndex1(indexRef) = spacing*(index-1)+1;
%     pointIndex2(indexRef) = spacing*(index-1)+1 + numPointsPerTurn;       % Off-by-one-turn
% 
%     indexRef;
% 
% end



spacing =resolution;
for index = 1:numTransverse
    %spacing = (Ratio_odd_integer+1);
    
    indexRef = index + numLinesAlongWire;
    x1Array(indexRef) = xArray(spacing*(index-1)+1); 
    x2Array(indexRef) = xArray(spacing*(index-1) + numPointsPerTurn+1);   % Off-by-one-turn 
    y1Array(indexRef) = yArray(spacing*(index-1)+1);
    y2Array(indexRef) = yArray(spacing*(index-1) + numPointsPerTurn+1);   % Off-by-one-turn
    z1Array(indexRef) = zArray(spacing*(index-1)+1);
    z2Array(indexRef) = zArray(spacing*(index-1) + numPointsPerTurn+1);   % Off-by-one-turn
    
    pointIndex1(indexRef) = spacing*(index-1)+1;
    pointIndex2(indexRef) = spacing*(index-1)+1 + numPointsPerTurn;       % Off-by-one-turn

    indexRef;

end
end 
