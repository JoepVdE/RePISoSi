% -------------------------------------------------------------------------
%  VIDEOplot.m  --  Build an MP4 of I(t) and Ic(t) per line element.
%  Part of RePISoSi - https://github.com/JoepVdE/RePISoSi  -  License: MIT
%  Author : J.L. Van den Eijnden, 2026
%           Original implementation by M. Mentink (Oct. 2022),
%           extended by J.L. Van den Eijnden and A. Vaskuri (Oct. 2023).
%
%  Required workspace variables:
%      IArrayhistory, IcArrayhistory, IcArray, tArrayhistory.
%
%  Output: Iarray_Isc_5turn240mmallon.mp4 in the current working directory.
%
%  TODO: parameterize the output filename / frame rate / y-axis limits.
% -------------------------------------------------------------------------

VIDEO_FILENAME = 'Iarray_Isc_5turn240mmallon.mp4';   % TODO: parameterize
FRAME_RATE     = 5;
YLIM_CURRENT_A = [-4000, 8000];

figure(2);
set(gcf, 'Color', 'w');
set(gcf, 'Position', [3, 4, 1500, 1000]);

v = VideoWriter(VIDEO_FILENAME, 'MPEG-4');
v.FrameRate = FRAME_RATE;
open(v);

cleanupV = onCleanup(@() close(v));   % ensure handle is closed even on error

for k = 1:size(IArrayhistory, 2)
    titleString = sprintf('t = %s s', num2str(tArrayhistory(k), 3));

    plot(IArrayhistory(:, k));
    hold on;
    plot(IcArrayhistory(:, k));
    hold off;

    title(titleString);
    legend('I in line element', 'I_c of line element', 'Location', 'southwest');
    grid on;
    ylim(YLIM_CURRENT_A);
    xlim([1, length(IcArray)]);
    ylabel('Current [A]');
    xlabel('Line element number [-]');

    frame = getframe(gcf);
    writeVideo(v, frame);
end
