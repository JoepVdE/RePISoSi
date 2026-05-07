% -------------------------------------------------------------------------
%  Conductivity.m  --  Stand-alone Bloch-Grueneisen sweep used during
%                       material-model exploration. Plots resistivity vs.
%                       temperature for three Debye temperatures.
%  Part of RePISoSi - https://github.com/JoepVdE/RePISoSi  -  License: MIT
%  Author : J.L. Van den Eijnden, 2026
%           Original implementation by M. Mentink (Oct. 2022),
%           extended by J.L. Van den Eijnden and A. Vaskuri (Oct. 2023).
%
%  This is a research scratchpad, not part of the production simulation
%  path. It produces a single comparison figure.
%
%  HISTORY: an earlier revision contained a stray syntax-error line
%      @(x)A*(x^(-1))).^n*(integral(...
%  which was unreachable but prevented the file from being parsed. The
%  faulty line has been removed.
% -------------------------------------------------------------------------

clear; close all; clc;

T  = (1:1:100).';     % Temperatures to sweep (K)
n  = 5;                % Bloch-Grueneisen exponent
s  = n;

% Physical constants / aluminium parameters
kf       = 1;          % Fermi-plane radius (arbitrary units, normalised)
Ef       = 11.7;       % Fermi energy (eV)
vf       = 2.02e8;     % Fermi velocity, aluminium (cm/s in original derivation)
kb       = 1.38e-23;   % Boltzmann constant (J/K)
Omega    = 1;          % Volume of the elementary cell (normalised)
e_charge = 1.6602e-19; % Elementary charge (C)
M_atom   = 27;         % Atomic mass aluminium (u)
theta_eff = 1;         % Coupling parameter (treat as 1 here)
z_offset  = 0;         % Residual resistivity offset

A = (4*pi/3) * kf^2 * Ef^2 * Omega ...
    / (e_charge^2 * kb * theta_eff * vf^2 * M_atom);

% Sweep over three Debye-like reference temperatures
thetaR_list = [400, 350, 300];
rho_curves  = zeros(numel(T), numel(thetaR_list));

for j = 1:numel(thetaR_list)
    thetaR = thetaR_list(j);
    for i = 1:length(T)
        rho_curves(i, j) = A * (T(i)/thetaR).^n ...
            .* integral(@(x) ((x.^s .* exp(x)) ./ ((exp(x) - 1).^2)), ...
                        0, thetaR./T(i)) + z_offset;
    end
end

figure;
plot(T, rho_curves, 'LineWidth', 1.2);
xlabel('Temperature [K]');
ylabel('\rho [\Omega \cdot m, normalised]');
legend(arrayfun(@(thR) sprintf('\\theta_R = %d K', thR), thetaR_list, ...
                'UniformOutput', false), 'Location', 'best');
grid on;
set(gcf, 'Color', 'w');
