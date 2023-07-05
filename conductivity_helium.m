function k_helium = conductivity_helium(T)
%CONDUCTIVITY_HELIUM in [W/(m*K)] 
%input: T in [K] and pressure in bar
%   Detailed explanation goes here

He = importdata("HeConductivity_0.01MPa.csv");
k_helium = interp1(He.data(:,1),He.data(:,2),T);
end

%%%references: A/AS 7" Technical Note 1334 (revised) Thermophysical
%%%Properties of Helium-4 from 0.8 to 1500 K with Pressures to 2000 MPa
%%%Vincent D. Arp Robert D. McCarty Daniel G. Friend

%linear interpolation of NIST data at 0.01 MPa (0.1 bar). Could be improved
%by adding data for more pressures. 