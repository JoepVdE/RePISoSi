clear
%I = 4500;
% T = 4:300;
% heat = heat_capacity_al_alloy_5083_burnout(T);
% 
% rho4K = BlochGrun(4,0.24233,432.43717,0.0247015,5);
% rhotable = BlochGrun(T,0.24233,432.43717,0.0247015,5)./rho4K;
% 
% R4K = 37E-7;
% R4300K = R4K.*rhotable;

% fun = @(T,t) (R/Cp).*I.^2;
% 
% q = integral(@(t) fun(T,))
% 
% 
% 
% tspan = [0 100];
% T0 = 4.2;
% [t,y] = ode45(@(t,T) 2*t, tspan, y0);

% syms Tt(t)
% ode = diff(Tt) == R()
% 
% ode = diff(y)+4*y == exp(-t);
% cond = y(0) == 1;
% ySol(t) = dsolve(ode,cond)
% integratedheatcapacity = trapz(heat);
% I = 4500;
% P = 4.5E-6*I^2;
% timeminutes= integratedheatcapacity/(P*60)
% plot(T,heat)
t_end = 1050;
steps = 5000;
Time = linspace(0,t_end,steps);
Initial_Conditions = [4.2];
[t,y] = ode45(@(t,y) Function_ODE1(t,y),Time,Initial_Conditions);


plot(t/60,y)
grid on
xlabel('t [min]')
ylabel('T [K]')
set(gcf,'color','w');
xlim([-10/60, t_end/60])

%% Function
function[dT_per_dt] = Function_ODE1(t,y)
    %weight = 0.5; %kg
    rho4K = BlochGrun(4,0.24233,432.43717,0.0247015,5);
    %R4K = 37E-7; %[Ohm] two times slots
    R4K = 9.8326E-6; %[Ohm] 
    I = 4500;
    Rf = R4K*BlochGrun(y,0.24233,432.43717,0.0247015,5)./rho4K;
    Cpf = heat_capacity_al_alloy_5083_burnout(y);
    dT_per_dt = (Rf/Cpf).*I.^2;
end

