T =  [1:1:100];
s = 5;              %power
n = s;

kf = 1;             %Radius of fermi plane
Ef = 11.7;          %[eV] Fermi energy
vf = 2.02E8;        %Fermi velocity Aluminum
kb = 1.38E-23;      %Boltzmann constant
Omega = 1;          % volume of elementarzelle
e = 1.6602E-19;     %elementary charge
M = 27;             %[u] Mass atom
thetaR = 400;       %[K] Debye temp
theta = 1;          % ????????????
z = 0;              % Start temp
A = (4*pi/3)*kf^2*Ef^2*Omega/(e^2*kb*theta*vf^2*M);

rhovalsout=zeros(size(T,2),size(T,1));
for i=1:length(T)
 rhovalsout(i)=A*(T(i)/thetaR).^n.*integral(@(x)((x.^s.*exp(x))./ ((exp(x)-1).^2)),0,thetaR./T(i))+z;
end

thetaR = 350;
rhovalsout2=zeros(size(T,2),size(T,1));
for i=1:length(T)
 rhovalsout2(i)=A*(T(i)/thetaR).^n.*integral(@(x)((x.^s.*exp(x))./ ((exp(x)-1).^2)),0,thetaR./T(i))+z;
end
thetaR = 300
rhovalsout3=zeros(size(T,2),size(T,1));
for i=1:length(T)
 rhovalsout3(i)=A*(T(i)/thetaR).^n.*integral(@(x)((x.^s.*exp(x))./ ((exp(x)-1).^2)),0,thetaR./T(i))+z;
end

%x = fminsearch(fun,x0)
%rho = factorial(n)*(zeta(n)-exp(-x)/(1-exp(-x))+ log(1-exp(-x))-symsum((x^j/(factorial(j))*polylog(n-j,exp(-x))),j,0,n-2))

plot(T,rhovalsout(:))
thetaR = 370;
for i=1:length(T)
@(x)A*(x^(-1))).^n*(integral(((x.^n.*exp(x))./ ((exp(x)-1).^2)),0,thetaR./T(i))+z
end


hold on
plot(T,rhovalsout2(:))
plot(T,rhovalsout3(:))
xlabel('Temperature [K]')
ylabel(' \rho [\Omega]')