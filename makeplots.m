% -------------------------------------------------------------------------
%  makeplots.m  --  Final summary plots for a completed run: centre B(t),
%                    coil E-field(t), and current-lead temperatures(t).
%  Part of RePISoSi - https://github.com/JoepVdE/RePISoSi  -  License: MIT
%  Author : J.L. Van den Eijnden, 2026
%           Original implementation by M. Mentink (Oct. 2022),
%           extended by J.L. Van den Eijnden and A. Vaskuri (Oct. 2023).
%
%  Required workspace variables (typically present after main_5turnMM.m
%  finishes or after loading a timestamped .mat dump):
%      tArrayhistory, centerBhistory, Pheatersave, Efield, Ecoil,
%      VGroundVectorhistory, VGroundVector,
%      temperatureRingTophistory, temperatureRingBottomhistory.
%
%  TODO: parameterize the hardcoded RSol = 0.12 m in the E-field tap
%        normalisation; replace with a workspace variable when available.
% -------------------------------------------------------------------------

figure(2);
set(gcf, 'Color', 'w');

% Voltage-tap probe indices used as a coil-averaged E-field estimate.
idxLow  = round(length(VGroundVector) * 1/10);   %#ok<NASGU>
idxHigh = round(length(VGroundVector) * 9/10);   %#ok<NASGU>

% Coil-averaged E-field (uV/m) per recorded iteration.
nFrames = size(VGroundVectorhistory, 2) - 1;
Ecoil   = zeros(1, nFrames);
nV      = length(VGroundVector);
for k = 1:nFrames
    Ecoil(k) = 1e6 * (VGroundVectorhistory(round(nV*1/10), k) ...
                    - VGroundVectorhistory(round(nV*9/10), k));
end

figure(7);
plot(tArrayhistory(1:length(Ecoil)), Ecoil);
xlabel('time [s]');
ylabel('E [\muV/m]');
grid on;
xlim([1, 400]);

% Centre-axis B(t) plus heater-power overlay.
figure(5);
plot(tArrayhistory(1:length(centerBhistory)), centerBhistory);
hold on;
plot(tArrayhistory(1:length(centerBhistory)), Pheatersave);
hold off;
xlabel('time [s]');
ylabel('B [T]');
grid on;
xlim([1, 400]);
ylim([0, 0.12]);
legend('B_{center}', 'P_{heater}', 'Location', 'best');

% Top/bottom current-lead temperatures.
figure(6);
plot(tArrayhistory(1:length(temperatureRingBottomhistory)), temperatureRingBottomhistory);
hold on;
plot(tArrayhistory(1:length(temperatureRingBottomhistory)), temperatureRingTophistory);
hold off;
grid on;
xlabel('time [s]');
ylabel('T [K]');
legend('T_{bottom lead}', 'T_{top lead}', 'Location', 'best');
xlim([1, 400]);
