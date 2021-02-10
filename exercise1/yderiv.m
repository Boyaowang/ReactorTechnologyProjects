function dydz=yderiv(z,y)
global rhob Ta R U dt dp molmass_so2 molmass_so3 ...
       molmass_o2 molmass_N2;
   
ptot = y(1);
T = y(2);
us = y(3);
pso2 = y(4);
po2 = y(5);
pso3 = y(6);
%pn2 = y(7); 
pn2 = ptot - pso2 - po2 -pso3;

% equation parameters

Mm_coeff = Mm(pso2,pso3,po2,pn2);

rhog=ptot*Mm_coeff/(R*T);

f_coeff = f(Rep(rhog,us));

Cpg = rho(pso2,molmass_so2,T) * CpSO2(T) ...
     + rho(po2,molmass_o2,T) * CpO2(T) ...
     + rho(pso3,molmass_so3,T) * CpSO3(T) ...
     + rho(pn2,molmass_N2,T) * CpN2(T);
enthalpy = deltaHR(T);

k_coeff = k(T);

rso2 = r(k_coeff,po2,pso2,pso3,Kp(T),"so2");
rso3 = r(k_coeff,po2,pso2,pso3,Kp(T),"so3");
ro2 = r(k_coeff,po2,pso2,pso3,Kp(T),"o2");

dTdz=1.0/(us*rhog*Cpg)*(-enthalpy*rhob*rso2-4.0*U/dt*(T-Ta));
dptdz=-f_coeff*rhog*us^2/dp; 
dusdz=-us/ptot*(dptdz-ptot/T*dTdz);%-us*Mm_coeff/(rhog*R) * (1/T*dptdz - ptot/(T^2)*dTdz);%
dpso2dz=-pso2/us*dusdz + pso2/T*dTdz - rhob*R*T*rso2/us;
dpo2dz=-po2/us*dusdz + po2/T*dTdz - rhob*R*T*ro2/us;
dpso3dz=-pso3/us*dusdz + pso3/T*dTdz - rhob*R*T*rso3/us;
%dpn2dz=-pn2/us*dusdz + pn2/T*dTdz;

dydz=[dptdz;dTdz;dusdz;dpso2dz;dpo2dz;dpso3dz];

end


%% heat capacity
function Cp = CpSO2(T)
Cp = 30.178+ T*42.452*10^(-3) - T*T*18.218*10^(-6);
end

function Cp = CpO2(T)
Cp = 23.995+ T*17.507*10^(-3) - T*T*6.628*10^(-6);
end

function Cp = CpSO3(T)
Cp = 35.634+ T*71.722*10^(-3) - T*T*31.539*10^(-6);
end

function Cp = CpN2(T)
Cp = 26.159+ T*6.615*10^(-3) - T*T*2.889*10^(-7);
end

%% Equilibrium constant
function Kp_ = Kp(T)
global R;
Kp_ = 3.142*10^(-3)*exp(98325/(R*T)-11.24); %[Pa^(-0.5)]
end

%% Reaction rare constant
function kRate = k(T)
kRate = 9.8692*10^(-3)*exp(-97782/T-110.1*log(T)+848.1); %[mol(SO2)/kg(cat)*Pa]
end

%% Reaction rate
function rrate = r(k,po2,pso2,pso3,Kp,name)
if(name == "so2")
        rrate = k* sqrt(abs(pso2/pso3)) * ( po2 - (pso3/(Kp * pso2))^2 );
end
if(name == "so3")
        rrate = -k* sqrt(abs(pso2/pso3)) * ( po2 - (pso3/(Kp * pso2))^2 );
end
if(name == "o2")
        rrate = 0.5*k* sqrt(abs(pso2/pso3)) * ( po2 - (pso3/(Kp * pso2))^2 );
end
end

%% Reynolds number
function Rey = Rep(rhog,us)
global dp mu;
Rey = rhog*us*dp/mu; %Calculating the friction factor;
end

%% Friction coefficient
function fri = f(Reynolds)
global epsilon;
fri = ((1-epsilon)/epsilon^3) * (1.75+ 4.2*Reynolds^(5/6)*(1-epsilon)/Reynolds);
end

%% reaction energy/enthalpy
function HR = deltaHR(T)
global Tr HR_Tr;
HR = HR_Tr + (-6.5415)* (T-Tr) + 0.5 * 0.02057 * (T^2 - Tr^2) + (-1.0011e-5)/3 * (T^3 - Tr^3);
end

%% Averaged molar mass
function Mavg = Mm(pso2,pso3,po2,pn2)
global molmass_so2 molmass_so3 molmass_o2 molmass_N2;
Mavg = (pso2 * molmass_so2 + pso3 * molmass_so3 ...
      + po2 * molmass_o2 + pn2 * molmass_N2) ...
      /(pso2 + pso3 + po2 + pn2);
end

%% phase density
function rhoi = rho(p,molmass,T)
global R;
rhoi = p * molmass / (R * T);
end
