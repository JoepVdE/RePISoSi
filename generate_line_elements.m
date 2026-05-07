% -------------------------------------------------------------------------
%  generate_line_elements.m  --  Generate the 3-D points of a uniform
%                                helical solenoid (cylindrical -> Cartesian).
%  Part of RePISoSi - https://github.com/JoepVdE/RePISoSi  -  License: MIT
%  Author : J.L. Van den Eijnden, 2026
%           Original implementation by M. Mentink (Oct. 2022),
%           extended by J.L. Van den Eijnden and A. Vaskuri (Oct. 2023).
%
%  Inputs : numPoints, numPointsPerTurn, RSolEff (effective radius, m),
%           numWindings, wCond (conductor width = pitch per turn, m).
%  Outputs: xArray, yArray, zArray  -- column vectors of Cartesian coords.
% -------------------------------------------------------------------------
function [xArray, yArray, zArray] = generate_line_elements(numPoints, ...
        numPointsPerTurn, RSolEff, numWindings, wCond)
    xArray = zeros(numPoints, 1);
    yArray = zeros(numPoints, 1);
    zArray = zeros(numPoints, 1);

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
