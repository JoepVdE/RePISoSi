% -------------------------------------------------------------------------
%  plotIcTc.m  --  Two-panel post-processing plot: |B|, Ic, |I| per element
%                  and the time evolution of E-field and helium-cooling load.
%  Part of RePISoSi - https://github.com/JoepVdE/RePISoSi  -  License: MIT
%  Author : J.L. Van den Eijnden, 2026
%           Original implementation by M. Mentink (Oct. 2022),
%           extended by J.L. Van den Eijnden and A. Vaskuri (Oct. 2023).
%
%  Run this AFTER a successful main_5turnMM.m simulation, or after loading
%  one of the timestamped .mat workspace dumps. Required workspace
%  variables:
%      BArray, IcArray, IArray, tArrayhistory, VGroundVectorhistory,
%      Qheliumcoolingarrayhistory.
%
%  TODO: parameterize the hardcoded RSol = 0.12 in the E-field denominator
%        and the iteration-count cap of 278 (currently used to size the
%        Qhelium time series).
% -------------------------------------------------------------------------

figure(2);
plot(BArray);
xlabel('Line element number [-]');
ylabel('|B| at line element [T]');
grid on;
set(gcf, 'Color', 'w');

figure(3);
plot(IcArray / 1000);
xlabel('Line element number [-]');
ylabel('I_c at line element [kA]');
grid on;
set(gcf, 'Color', 'w');

figure(4);
plot(IArray / 1000);
xlabel('Line element number [-]');
ylabel('I in line element [kA]');
grid on;
set(gcf, 'Color', 'w');

% E-field over the 5-turn helix (probe = first voltage tap).
% TODO: parameterize - the (2*pi*0.12) below assumes RSol = 0.12 m.
figure(5);
plot(tArrayhistory, 1e6 * VGroundVectorhistory(1, :) / (2 * pi * 0.12));
xlabel('time [s]');
ylabel('E over 5 turns [\muV/m]');
xlim([0, 50]);
grid on;

% Total helium-cooling load history.
nIter   = size(Qheliumcoolingarrayhistory, 3);
Qhelium = zeros(1, nIter);
for i = 1:nIter
    Qhelium(i) = sum(sum(Qheliumcoolingarrayhistory(:, :, i)));
end

figure(6);
plot(tArrayhistory(1:nIter), Qhelium);
xlabel('time [s]');
ylabel('Q helium [W]');
grid on;
