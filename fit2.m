clear 
M = importdata('Alres.csv');

Tf=400;
%T2=M.data(:,1);
%rho=M.data(:,2)*1E-4;
T=M.data(:,3);
rho2=M.data(:,4)*1E-8;


ft = fittype('BlochGrun(x,A,thetaR,z,n)')
options = fitoptions(ft);
options.Startpoint = [2E-7, 430,1E-12,5]
options.Lower = [0,10,0,4.5]
options.Upper = [0.1,1000,1,5.5]
%T=x;
x=T;
f = fit(x, rho2, ft, options )
plot(f,x,rho2) 
%Upper[A]=1
grid on
xlabel('T[K]')
ylabel('\rho [\Omega m]')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')