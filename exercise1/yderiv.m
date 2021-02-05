function dydz=yderiv(z,y)
global molmass rhob pb0 Cp enthalpy U dt Tr R;
pso2 = y(1);
po2 = y(2);
pso3 = y(3);
T = y(4);
ptot = y(5);
us = y(6);
f=1;
dp=0.1;
k=1.0e-10*exp(19.837-13636.0/T)/3600.0; 
r=k*pb0*pso2;
rso2 = r;
ro2 = 0.5 * r;
rso3 = -r;
rhog=ptot*molmass/(R*T);
dptdz=-f*rhog*us^2/dp;
dTdz=1.0/(us*rhog*Cp)*(enthalpy*rhob*r-4.0*U/dt*(T-Tr));
dusdz=-us*molmass/(rhog*R)*(1/T*dptdz-ptot/(T^2)*dTdz);
dpso2dz=-1.0/us*(molmass*ptot*rhob/rhog*r)+pso2/T*dTdz+pso2/us*dusdz;
dpo2dz=-1.0/us*(molmass*ptot*rhob/rhog*r)+po2/T*dTdz+po2/us*dusdz;
dpso3dz=-1.0/us*(molmass*ptot*rhob/rhog*r)+pso3/T*dTdz+pso3/us*dusdz;

dydz=[dptdz;dTdz;dusdz;dpso2dz];
end


