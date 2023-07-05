
function [k] = thermal_conductivity_Al10SiMg(T)

    
%Wiedemann–Franz law
L = 2.44E-8; %[V^2K^-2]
k = 1E6.*T.*L./BlochGrun(T,0.24233,432.43717,0.0247015,5);

    
end
