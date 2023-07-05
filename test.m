clear
% 
% syms m
% n=5;
% 
% Tf  = 400; %[K] fermi temperature
% T = 1:1:370;    %[K] temp array
% x = Tf./T;
% 
% 
% %x = 0.01:0.001:1;
% 
% integral = factorial(n)*(zeta(n)-exp(-x)/(1-exp(-x))+ log(1-exp(-x))-symsum(((x.^m./factorial(m)).*polylog(n-m,exp(-x))),m,0,n-2));
% rho = (x.^(-n)).*integral;
% double(rho)
% plot(x,rho)
T = 1:100;
thetaR = 370;
n = 5;
for i=1:length(T)
rho = @(x)(A*(x^(-1))).^n*(integral(((x.^n.*exp(x))./ ((exp(x)-1).^2)),0,thetaR./T(i))+z)
end
plot(T,rho(T))