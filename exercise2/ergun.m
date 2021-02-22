% ergun
% This function computes the 1D pressure gradient in the system by using
% superficial velocity, density and Reynolds number.
% Input:
% RHO       [=] kg/m^3          Gas density
% u         [=] m/s             Axial velocity
% Re        [=] -               Reynolds number
% r         [=] m               Radial positioning vector
%
% Output
% dpdz      [=] Pa/m            Pressure gradient

function dpdz = ergun(RHO,u,Re,r)

global EPS Dp RP

%Area averaged density, velocity and Reynolds number for the Ergun equation
SUMRHO = 0;
SUMu = 0;
SUMRe = 0;
SUMRHO = SUMRHO + RHO(1)*r(1)^2;
SUMu = SUMu + u(1)*r(1)^2;
SUMRe = SUMRe + Re(1)*r(1)^2;
for i=2:RP
    SUMRHO = SUMRHO + RHO(i)*(r(i)^2-r(i-1)^2);
    SUMu = SUMu + u(i)*(r(i)^2-r(i-1)^2);
    SUMRe = SUMRe + Re(i)*(r(i)^2-r(i-1)^2);
end
RHOav = SUMRHO/r(RP)^2;
uav = SUMu/r(RP)^2;
Reav = SUMRe/r(RP)^2;

%Friction factor
a = 1.75;
b = 4.2*Reav.^(5/6);
f = (1-EPS)/EPS^3*(a+b*(1-EPS)/Reav);

%The Ergun equation
dpdz = -f*RHOav*uav^2/Dp;