clc;
close all;
clear all;

global rhob Ta R U dt dp mu epsilon molmass_so2 molmass_so3 ...
       molmass_o2 molmass_N2 Tr HR_Tr;

rhob = 541.42; %[kg/m^3]
Ta = 702.6; %[K]
R = 8.3145; %[J/mole*K]
U = 56.783; %[J/(m^2*s*K)]
dt = 2*0.0353; %[m]
dp = 0.004572;  %[m]
mu = 3.7204e-5; %[kg/(m*s)]
epsilon = 0.45;
molmass_so2 = (15.999*2+32.066) / 1000; %[kg/mole] 
molmass_so3 = (15.999*3+32.066) / 1000; %[kg/mole]
molmass_o2 = (15.999*2) / 1000; %[kg/mole]
molmass_N2 = (14.007*2) / 1000; %[kg/mole]
Tr = 699.8; %[K]
HR_Tr = -98787.5; %[J/(molSO2)]

% initial conditions
nt0 = 54.8214; %[mol/(m2 s)]
T0 = 777.78; %[K]
p0 = 202650; %[pa]

us_0 = nt0*R*T0/p0; %[m/s] initial superfacial velocity
pso2_0 = 22291.5; %[pa]
pso3_0 = eps; %[pa]
po2_0 = 20265; %[pa]
%pn2_0 = p0-pso2_0-pso3_0-po2_0; %[pa]

zstart = 0; %[m]
zend = 6.096; %[m]

%integration span
zspan=[zstart zend];

%initial conditions
y0=[p0 T0 us_0 pso2_0 po2_0 pso3_0];

[z,y]=ode15s(@yderiv,zspan,y0);

pn2Plot = y(:,1) - y(:,4) -y(:,5) -y(:,6);

%plot the result

m = 3;
n = 2;
nr = 1;

subplot(m,n,nr);
plot(z,y(:,1))
title('total pressure profile')
xlabel('z [m]') 
ylabel('Ptot [Pa]')

subplot(m,n,2);
plot(z,y(:,2))
title('Temperature profile')
xlabel('z [m]') 
ylabel('T [K]')

subplot(m,n,3);
plot(z,y(:,3))
title('superfacial velocity profile')
xlabel('z [m]') 
ylabel('us [m/s]')

subplot(m,n,4);
plot(z,y(:,4))
title('SO2 pressure profile')
xlabel('z [m]') 
ylabel('PSO2 [pa]')

subplot(m,n,5);
plot(z,y(:,5))
title('O2 pressure profile')
xlabel('z [m]') 
ylabel('PO2 [pa]')

subplot(m,n,6);
plot(z,y(:,6))
title('SO3 pressure profile')
xlabel('z [m]') 
ylabel('PSO3 [pa]')

% subplot(m,n,7);
% plot(z,pn2Plot)
% title('N2 pressure profile')
% xlabel('z [m]') 
% ylabel('PN2 [pa]')

figure()
plot(z,y(:,1))
hold on
plot(z,y(:,4))
hold on
plot(z,y(:,5))
hold on
plot(z,y(:,6))
hold on
plot(z,pn2Plot)

legend({'Ptot','PSO2','PO2', 'PSO3','PN2'},'Location','northeast')









