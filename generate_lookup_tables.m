tic
clear
close all
fprintf('clear \n')
toc
Tinitialise = 0.1:0.1:400;



L = 2.44E-8; %[V^2K^-2]
ktable = 1E6.*Tinitialise.*L./BlochGrun(Tinitialise,0.2594,484,0.026,5);
rhotable = 1E6.*BlochGrun(Tinitialise,0.2594,484,0.026,5);

rng(1); %get repeatable set of random numbers
temparray = 5*ones(600,5);
fprintf('make tables \n')
toc 
for j = 1 : size(temparray,2)
for i = 1 : size(temparray,1)
    temparray(i,j)= temparray(i,j)+4*(rand-0.5);
end
end
fprintf('random numbers \n')
toc 
rounded = round(10*temparray(:,:));
kmatrix = zeros(size(temparray));
kmatrix = ktable(rounded(:,:));
rhomatrix = zeros(size(temparray));
rhomatrix = rhotable(rounded(:,:));
% 
% figure(1)
% plot(Tinitialise,ktable)
% grid on
% 
% figure(2)
% plot(Tinitialise,rhotable)
% grid on
fprintf('get values \n')
toc
