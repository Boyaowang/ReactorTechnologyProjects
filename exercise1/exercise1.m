clc;
clear all;

global rhob Ta R U dt dp mu epsilon;

nt0 = 54.8214; %[mol/(m2 s)]
supvel_0 = 2.2e-3; %[m/s] initial superfacial velocity
molmass = 29.48; %[kg/kmole] 
ptot = 1.0e5; %[Pa] 
rhob = 541.42; %[kg/m^3]
Ta = 702.6; %[K]
pb0 = 0.211e5; %[Pa]
Cp = 0.992; %[kJ/kg*K] 
enthalpy = 1285409.0; %[kJ/kmole]
U = 56.783; %[J/m^2*s*K]
dt = 2*0.0353; %[m]
dp = 0.004572;  %[m]
mu = 

R = 8.3145; %[J/mole*K]
zstart = 0; %[m]
zend = 3; %[m]
pA0 = 0.015e5; %[Pa]
T0 = 625; %[K]
ptot0 = 10000;
us0 = 2.2e-3;

%integration span
zspan=[zstart zend];

%initial conditions
y0=[pA0 T0 ptot0 us0];

[z,y]=ode15s(@yderiv,zspan,y0);

%plot the result

m = 4;
n = 1;
nr = 1;


subplot(m,n,nr);
plot(z,y(:,1))
title('Partial pressure profile')
xlabel('z [m]') 
ylabel('p [Pa]')

subplot(m,n,2);
plot(z,y(:,2))
title('Temperature profile')
xlabel('z [m]') 
ylabel('T [K]')


subplot(m,n,3);
plot(z,y(:,3))
title('total pressure profile')
xlabel('z [m]') 
ylabel('Ptot [pa]')


subplot(m,n,4);
plot(z,y(:,4))
title('superfacial velocity profile')
xlabel('z [m]') 
ylabel('us [m/s]')



