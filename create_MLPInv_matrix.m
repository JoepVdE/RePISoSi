function MLPInv = create_MLPInv_matrix(MArray, numLines, numPoints, pointIndex1, pointIndex2)
    % Create matrix that keeps track of equations for lines and points
    MLP = zeros(numLines + numPoints - 1, numLines + numPoints - 1);
    MLP(1:numLines, 1:numLines) = MArray;

    % Build the matrix equations for the lines
    % The inductive voltage contribution is considered
    % Subsequently, the voltage drop is considered: For example between points
    % 2 and 3, the voltage of point 3 is added and that of point 2 is
    % substracted
    % Note that point 1 is always at 0 V (i.e. it is biased), and therefore not
    % present in there equations.

    for index = 1:numLines
        pointPlusIndex = pointIndex1(index) + numLines;
        pointMinIndex = pointIndex2(index) + numLines;

        if(pointPlusIndex < (numLines + numPoints))
            MLP(index, pointPlusIndex) = -1;
        end
        if(pointMinIndex < (numLines + numPoints))
            MLP(index, pointMinIndex) = 1;
        end

    end

    % Build the equations for the points
    for index = 1:(numPoints - 1)
        indexPoint = index + numLines;

        for index2 = 1:numLines
            % Checks per point whether current is flowing out or into a point
            if(pointIndex1(index2) == index)
                MLP(indexPoint, index2) = -1;
            end
            if(pointIndex2(index2) == index)
                MLP(indexPoint,index2) = 1;
            end
        end    
    end    

    % ***** Matrix inversion ***** 
    % NOTE: 
    % Y = inv(X) computes the inverse of square matrix X.
    % X^(-1) is equivalent to inv(X).
    % x = A\b is computed differently than x = inv(A)*b and is recommended for solving systems of linear equations.
    % See: https://www.mathworks.com/help/matlab/ref/inv.htm
    disp("Matrix inversion");
    f = waitbar(0, 'Matrix inversion');
    MLPInv = MLP^-1; 
    close(f);
end
