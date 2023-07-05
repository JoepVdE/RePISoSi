function [xArray, yArray, zArray] = generate_line_elements(numPoints, numPointsPerTurn, RSolEff, numWindings, wCond) 
    % Generate line elements for solenoid 
    xArray = zeros(numPoints, 1); % Initialize points in x-axis
    yArray = zeros(numPoints, 1); % Initialize points in y-axis
    zArray = zeros(numPoints, 1); % Initialize points in z-axis 

    len = numWindings*wCond; 
    dzdp = len/(numPointsPerTurn*numWindings); % dz along z-axis of the solenoid 

    pause(0.1);
    tic

    f = waitbar(0, 'Generating line elements');
    updateCounter = floor(numPoints/100);
    for index = 1:numPoints 
        % Conversion from cylindrical to cartesian coordinate system 
        x = RSolEff*cos(2*pi*(index - 1)/(numPointsPerTurn));
        y = RSolEff*sin(2*pi*(index - 1)/(numPointsPerTurn));
        z = -len/2 + dzdp*(index - 1);   

        % Generate line elements
        xArray(index) = x; % Along x-axis
        yArray(index) = y; % Along y-axis 
        zArray(index) = z; % Along z-axis    

        if(mod(index, updateCounter) == 0)
            waitbar(index/numPoints, f, 'Generating line elements');
            pause(0.1);
        end    
    end
    close(f);
end 
