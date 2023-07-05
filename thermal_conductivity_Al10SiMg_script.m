
Tarray = 4:0.1:300;

k = thermal_conductivity_Al10SiMg(Tarray);

figure(4)
plot(Tarray,k,LineWidth=2)
grid on
xlabel('{\it T} [K]')
ylabel('{\it k} [ W/(m\cdotK)]')
hold on
plot(Tarray,thermal_conductivity_al_alloy_5083(Tarray),LineWidth=2)
hold off

legend('Al10SiMg','Al5083')
