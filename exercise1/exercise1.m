clc;
clear all;

global rhob Ta R U dt dp mu epsilon molmass_so2 molmass_so3 ...
       molmass_o2 molmass_N2;

rhob = 541.42; %[kg/m^3]
Ta = 702.6; %[K]
R = 8.3145; %[J/mole*K]
U = 56.783; %[J/(m^2*s*K)]
dt = 2*0.0353; %[m]
dp = 0.004572;  %[m]
mu = 3.7204e-5; %[kg/(m*s)]
epsilon = 0.45;
molmass_so2 = 15.999*2+32.066; %[kg/kmole] 
molmass_so3 = 15.999*3+32.066; %[kg/kmole]
molmass_o2 = 15.999*2; %[kg/kmole]
molmass_N2 = 14.007*2; %[kg/kmole]

% initial conditions
nt0 = 54.8214; %[mol/(m2 s)]
T0 = 777.78; %[K]
p0 = 202650; %[pa]
us_0 = nt0*R*T0/p0; %[m/s] initial superfacial velocity
pso2

ptot = 1.0e5; %[Pa] 
pb0 = 0.211e5; %[Pa]
Cp = 0.992; %[kJ/kg*K] 
enthalpy = 1285409.0; %[kJ/kmole]


zstart = 0; %[m]
zend = 3; %[m]
pA0 = 0.015e5; %[Pa]
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



