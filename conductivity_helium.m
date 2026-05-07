% -------------------------------------------------------------------------
%  conductivity_helium.m  --  Linear interpolation of Helium-4 thermal
%                              conductivity at p = 0.01 MPa (0.1 bar).
%  Part of RePISoSi - https://github.com/JoepVdE/RePISoSi  -  License: MIT
%  Author : J.L. Van den Eijnden, 2026
%           Original implementation by M. Mentink (Oct. 2022),
%           extended by J.L. Van den Eijnden and A. Vaskuri (Oct. 2023).
%
%  Input : T (K)        - any shape, must lie within the table range.
%  Output: k_helium (W/(m*K))
%
%  Reference: Arp, McCarty, Friend, NBS Technical Note 1334 (rev.),
%             Thermophysical Properties of Helium-4 from 0.8 to 1500 K
%             with Pressures to 2000 MPa.
%
%  TODO: parameterize - currently hardcoded to 0.01 MPa data; extend to
%        multi-pressure tables if helium pressure becomes a free parameter.
% -------------------------------------------------------------------------
function k_helium = conductivity_helium(T)
    He = importdata("HeConductivity_0.01MPa.csv");   % TODO: parameterize file path
    k_helium = interp1(He.data(:,1), He.data(:,2), T);
end