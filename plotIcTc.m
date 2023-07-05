
clear all; close all; clc; 

set(0, 'defaultTextFontName', 'times', 'defaultTextFontSize', 14);
set(0, 'defaultAxesFontName', 'times', 'defaultAxesFontSize', 14);

T = (0:0.1:30)'; % Temperature 

Tc = 9; % Critical temperature 
Ic = 200*(Tc - T)./Tc; 
Ic = Ic.*(Ic > 0); % Critical current as a function of temperature 

figure(1);
set(gcf, 'color', 'w'); 
plot(T, Ic); 
xlabel('Temperature / K'); 
ylabel('Critical current / A'); 

