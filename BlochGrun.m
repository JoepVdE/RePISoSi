% -------------------------------------------------------------------------
%  BlochGrun.m  --  Bloch-Grueneisen resistivity formula (vectorised in T).
%  Part of RePISoSi - https://github.com/JoepVdE/RePISoSi  -  License: MIT
%  Author : J.L. Van den Eijnden, 2026
%           Original implementation by M. Mentink (Oct. 2022),
%           extended by J.L. Van den Eijnden and A. Vaskuri (Oct. 2023).
%
%  rho = BlochGrun(T, A, thetaR, z, n) returns the temperature-dependent
%  resistivity for an array of temperatures T (K). Parameters:
%       A      : prefactor (Ohm*m or scaled, see calling site)
%       thetaR : Debye-like reference temperature (K)
%       z      : residual resistivity offset
%       n      : exponent (typically 5 for clean metals)
%  Inputs T < TLowLimit = 2 K are clamped to TLowLimit to avoid the
%  integral becoming numerically singular.
% -------------------------------------------------------------------------
function rho = BlochGrun(x, A, thetaR, z, n)
    TLowLimit = 2;
    rhovalsout = zeros(size(x,2), size(x,1));




    % Clamp very low temperatures (added by M. Mentink, 28/07/2023) to avoid
    % numerical warnings from the integral kernel near T = 0.
    x(x < TLowLimit) = TLowLimit;

    for i = 1:length(x)
        rhovalsout(i) = A .* (x(i)./thetaR).^n ...
            .* integral(@(u) ((u.^n .* exp(u)) ./ ((exp(u) - 1).^2)), ...
                        0, thetaR./x(i)) + z;
    end
    rho = rhovalsout';
end
