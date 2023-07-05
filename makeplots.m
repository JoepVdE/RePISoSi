figure(2)
plot(BArray)
xlabel('Line element number [-]')
ylabel('|B| at line element [T]')
grid on
set(gcf,"Color",'w')
figure(3)
plot(IcArray/1000)
xlabel('Line element number [-]')
ylabel('I_c at line element [kA]')
grid on
set(gcf,"Color",'w')
figure(4)
plot(IArray/1000)
xlabel('Line element number [-]')
ylabel('I in line element [kA]')
grid on

set(gcf,"Color",'w')


figure(5)
 plot(tArrayhistory,1e6*VGroundVectorhistory(1,:)/(2*pi*0.12))
 xlabel('time [s]')
ylabel('E over 5 turns [\muV/m]')
xlim([0 50])
grid on

for i = 1:278
Qhelium(i) = sum(sum(Qheliumcoolingarrayhistory(:,:,i)));
end
figure(6)
plot(tArrayhistory,Qhelium)
 xlabel('time [s]')
ylabel('Q helium [W]')
