% -------------------------------------------------------------------------
%  plotBET.m  --  Snapshot E-field probe over the 5-turn helix at a chosen
%                  iteration. Originally also produced a Vground.mp4 video
%                  but the recording loop is left disabled by default.
%  Part of RePISoSi - https://github.com/JoepVdE/RePISoSi  -  License: MIT
%  Author : J.L. Van den Eijnden, 2026
%           Original implementation by M. Mentink (Oct. 2022),
%           extended by J.L. Van den Eijnden and A. Vaskuri (Oct. 2023).
%
%  Required workspace variables:
%      VGroundVectorhistory, tArrayhistory.
%
%  Set MAKE_VIDEO = true to re-enable the per-frame video writer.
%  TODO: parameterize the hardcoded RSol = 0.12 m used in the normalisation.
% -------------------------------------------------------------------------

MAKE_VIDEO = false;   % set true to write Vground5.mp4
RSol_assumed = 0.12;  % m  -- TODO: parameterize / read from workspace

figure(2);
set(gcf, 'Color', 'w');

if MAKE_VIDEO
    v = VideoWriter('Vground5.mp4', 'MPEG-4');   % TODO: parameterize output path
    v.FrameRate = 1;
    open(v);

    for k = 1:size(VGroundVectorhistory, 2)
        if tArrayhistory(k) > 40
            break;
        end
        plot(1e6 .* VGroundVectorhistory(:, k) ./ (5 * 2*pi*RSol_assumed));
        title(sprintf('t = %s s', num2str(tArrayhistory(k), 3)));
        legend('E [\muV/m]', 'Location', 'southwest');
        grid on;
        ylabel('E [\muV/m]');
        xlabel('line element nr [-]');
        frame = getframe(gcf);
        writeVideo(v, frame);
    end
    close(v);
end

% Single snapshot (last available iteration) when video is disabled.
k = size(VGroundVectorhistory, 2);
plot(1e6 .* VGroundVectorhistory(:, k) ./ (5 * 2*pi*RSol_assumed));
title(sprintf('Snapshot at t = %s s', num2str(tArrayhistory(k), 3)));
xlabel('line element nr [-]');
ylabel('E [\muV/m]');
grid on;
