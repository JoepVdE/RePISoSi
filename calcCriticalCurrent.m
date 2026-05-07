% -------------------------------------------------------------------------
%  calcCriticalCurrent.m  --  Quick visualisation of the linear Ic(T) law.
%  Part of RePISoSi - https://github.com/JoepVdE/RePISoSi  -  License: MIT
%  Author : J.L. Van den Eijnden, 2026
%           Original implementation by M. Mentink (Oct. 2022),
%           extended by J.L. Van den Eijnden and A. Vaskuri (Oct. 2023).
%
%  Plots Ic(T) = Ic0 * (Tc - T) / Tc, clamped to non-negative values, for
%  the educational/illustrative parameters Ic0 = 200 A and Tc = 9 K.
%  This script is not part of the production solver path.
% -------------------------------------------------------------------------

clear; close all; clc;

set(0, 'defaultTextFontName', 'times', 'defaultTextFontSize', 14);
set(0, 'defaultAxesFontName', 'times', 'defaultAxesFontSize', 14);

T   = (0:0.1:30).';     % Temperature sweep (K)
Tc  = 9;                 % Critical temperature (K)
Ic0 = 200;               % Critical current near 0 K (A)

Ic = Ic0 * (Tc - T) ./ Tc;
Ic = Ic .* (Ic > 0);     % clamp to >= 0

figure(1);
set(gcf, 'color', 'w');
plot(T, Ic, 'LineWidth', 1.2);
xlabel('Temperature / K');
ylabel('Critical current / A');
grid on;
