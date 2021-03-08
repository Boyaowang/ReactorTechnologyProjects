clc;
close all;
clear all;
%%%%%%%%%%%%%%%%%%% Loading the constants %%%%%%%%%%%%%%%%
constant
eta = 0.001; % efficiency factor

%%%%%%%%%%%%%%%%% Sizing the r- and z-direction %%%%%%%%%%%%%
% r-direction:
r0 = 0; % Lower integration limit in r-direction
rn = RADIUSi; % Upper integration limit in r-direction
r = (r0:(rn-r0)/(RP-1):rn)';

%%%%%%%%%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%
wCH40 = FRACin(1) .*ones(RP,1);
wCO0 = FRACin(2) .*ones(RP,1);
wCO20 = FRACin(3) .*ones(RP,1);
wH20 = FRACin(4) .*ones(RP,1);
wH2O0 = FRACin(5) .*ones(RP,1);
wN2 = FRACin(6) .*ones(RP,1);
T0 = Tin*ones(RP,1);

%initial guess vectors:
init = [wCH40; wCO0; wCO20; wH20; wH2O0; T0];
%parameters vectors:
par = [r0 eta uin];
%%%%%%%%%%%%%%%%%%%%% non-linear system solver %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create ODE solver:
options = optimoptions('fsolve','Display','none','PlotFcn',@optimplotfirstorderopt);

func = @calc;
results = fsolve(@calc, init, options);



% %Define the result vector:
% n = RP;
% wCH4 = y(:,1:n);
% wCO = y(:,n+1:2*n) ;
% wCO2 = y(:,2*n+1:3*n);
% wH2 = y(:,3*n+1:4*n);
% wH2O =y(:,4*n+1:5*n);
% wN2 = ones(size(wCH4)) - wCH4 - wCO -wH2 -wH2O;
% T = y(:,5*RP+1:6*RP);
% %%%%%%%%%%%%%%%%%%%%%%% Ploting the results %%%%%%%%%%%%%%%%%%%%
% 
% %plot the tempt + mass fraction profiles:
% for i=0:5
% plottingFunc(i,r,z,y,RP);
% subplot(3,2,6);
% mesh(r,z,wN2)
% axis([0 RADIUSi 0 LENGTH])
% grid on
% xlabel('Radius [m]')
% ylabel('Z [m]')
% zlabel('N2 mass fraction [-]');
% end
% 
% 
% %Plot the pressure profile
% figure()
% plot(z,y(:,(Ncomp+1)*RP+1))
% title('total pressure profile')
% xlabel('z [m]') 
% ylabel('Ptot [Pa]')
% 
% % % Plotting the Temperature
% figure
% mesh(r,z,T)
% grid on
% xlabel('Radius [m]')
% ylabel('Z [m]')
% zlabel('Temperature [K]');
