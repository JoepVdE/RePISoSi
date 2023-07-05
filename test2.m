clear 
M = importdata('Alres.csv');

Tf=400;
T=M.data(:,1);
rho=M.data(:,2)*1E-4;
T2=M.data(:,3);
rho2=M.data(:,4)*1E-8;

A = 1;
Tf = 360;
n = 5;
x = [4:300];
rho0 = 0;
rho = A.*(x./Tf).^(n).*(120).*( zeta(n) - (exp(-(Tf./x))./(1-exp(-(Tf./x)))) + log(1-exp(-(Tf./x))) ...
    - polylog(5,exp(-(Tf./x))) ...
    - (Tf./x).*polylog(4,exp(-(Tf./x))) ...
    - ((Tf./x).^2)./2.*polylog(3,exp(-(Tf./x))) ... 
    - ((Tf./x).^3)./6.*polylog(2,exp(-(Tf./x)))) + rho0;


rho2 = A.*(x./Tf).^(n).*(120).*( zeta(n) - (exp(-(Tf./x))./(1-exp(-(Tf./x)))) + log(1-exp(-(Tf./x))) ...
    - (Tf./x).^3)./6.*polylog(5,exp(-(Tf./x))) ...
    - ((Tf./x).^2)./2.*polylog(4,exp(-(Tf./x))) ...
    - (Tf./x).*polylog(3,exp(-(Tf./x))) ... 
    - (polylog(2,exp(-(Tf./x)))) + rho0;


rho3 = A.*(x./Tf).^(n).*(120).*( zeta(n) + (exp(-(Tf./x))./(1-exp(-(Tf./x)))) + log(1-exp(-(Tf./x))) ...
    + (polylog(5,exp(-(Tf./x))) ...
    + ((Tf./x).*polylog(4,exp(-(Tf./x))) ...
    + ((Tf./x).^2)./2.*polylog(3,exp(-(Tf./x))) ... 
    + ((Tf./x).^3)./6.*polylog(2,exp(-(Tf./x)))))) + rho0;


plot(x,rho)
hold on
plot(x,rho2)
plot(x,rho3)
xlabel('T [K]')
ylabel('\rho [-]')
grid on
legend("1","2","3")
hold off