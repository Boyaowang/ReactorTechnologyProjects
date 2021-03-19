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
ndisk = RP; % Number of discretization points
r = (r0:(rn-r0)/(ndisk-1):rn)';

% z-direction:
zstart = 0; %[m]
zend = LENGTH; %[m]
zspan=[zstart zend]; %integration span

%%%%%%%%%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%
wCH40 = FRACin(1) .*ones(ndisk,1);
wCO0 = FRACin(2) .*ones(ndisk,1);
wCO20 = FRACin(3) .*ones(ndisk,1);
wH20 = FRACin(4) .*ones(ndisk,1);
wH2O0 = FRACin(5) .*ones(ndisk,1);
wN2 = FRACin(6) .*ones(ndisk,1);
T0 = Tin*ones(ndisk,1);
uz0 = uin*ones(ndisk,1);
pt0 = pin;
%initial condition vectors:
init = [wCH40; wCO0; wCO20; wH20; wH2O0; T0; uz0; pt0];
%parameters vectors:
par = [r0 eta uin];
%%%%%%%%%%%%%%%%%%%%% ODEs solver %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create Mass matrix:
M = eye((Ncomp+1)*RP+1);
for i=1:ndisk:(Ncomp)*ndisk+1
    M(i,i) = 0;
    M(i+RP-1,i+RP-1) =0;
end
% Create ODE solver:
options = odeset('Mass', M, 'Stats', 'on');
% Solving the system of equations
[z,y] = ode15s(@calcGas,zspan,init,options, par,r);

%Define the result vector:
n = RP;
wCH4 = y(:,1:n);
wCO = y(:,n+1:2*n) ;
wCO2 = y(:,2*n+1:3*n);
wH2 = y(:,3*n+1:4*n);
wH2O =y(:,4*n+1:5*n);
wN2 = ones(size(wCH4)) - wCH4 - wCO -wH2 -wH2O;
T = y(:,5*ndisk+1:6*ndisk);
%%%%%%%%%%%%%%%%%%%%%%% Ploting the results %%%%%%%%%%%%%%%%%%%%

%plot the tempt + mass fraction profiles:
for i=0:5
plottingFunc(i,r,z,y,ndisk);
subplot(3,2,6);
mesh(r,z,wN2)
axis([0 RADIUSi 0 LENGTH])
grid on
xlabel('Radius [m]')
ylabel('Z [m]')
zlabel('N2 mass fraction [-]');
end


%Plot the pressure profile
figure()
plot(z,y(:,(Ncomp+1)*ndisk+1))
title('total pressure profile')
xlabel('z [m]') 
ylabel('Ptot [Pa]')

% % Plotting the Temperature
figure
mesh(r,z,T)
grid on
xlabel('Radius [m]')
ylabel('Z [m]')
zlabel('Temperature [K]');
