% -------------------------------------------------------------------------
%  analyseaxialcurrent.m  --  Visualise axial-short currents on top of the
%                              helix, for a loaded simulation result.
%  Part of RePISoSi - https://github.com/JoepVdE/RePISoSi  -  License: MIT
%  Author : J.L. Van den Eijnden, 2026
%           Original implementation by M. Mentink (Oct. 2022),
%           extended by J.L. Van den Eijnden and A. Vaskuri (Oct. 2023).
%
%  Required workspace variables:
%      x1Array, x2Array, y1Array, y2Array, z1Array, z2Array, numPoints
%
%  Plots the helix points (small blue dots) and overlays the axial-short
%  segments (red-filled markers).
% -------------------------------------------------------------------------

figure(2);

% Pre-allocate plot arrays for axial-short segments.
nShortPoints = 2 * (length(x1Array) - numPoints + 1);
plotarrayx = zeros(nShortPoints, 1);
plotarrayy = zeros(nShortPoints, 1);
plotarrayz = zeros(nShortPoints, 1);

j = 1;
for i = numPoints:length(x1Array)
    plotarrayx(j)     = x1Array(i);
    plotarrayx(j + 1) = x2Array(i);
    plotarrayy(j)     = y1Array(i);
    plotarrayy(j + 1) = y2Array(i);
    plotarrayz(j)     = z1Array(i);
    plotarrayz(j + 1) = z2Array(i);
    j = j + 2;
end

plot3(x1Array, y1Array, z1Array, 'o', 'Color', 'b', 'MarkerSize', 1);
hold on;
plot3(plotarrayx, plotarrayy, plotarrayz, '-o', ...
    'Color', 'b', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
hold off;
grid on;
axis equal;
xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');
title('Helix mesh points (blue) and axial-short segments (red)');
