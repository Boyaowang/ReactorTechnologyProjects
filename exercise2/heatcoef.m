% heatcoef
% The function computes the heat transfer coefficient for radial transport of
% heat from the bed to the surroundings
% Input:
% Re        [=] -               Reynolds number
% T         [=] K               Temperature
% Y         [=] -               Mol fraction
% VIS       [=] kg/ms           Gas viscosity
% CPgas     [=] J/kgK           Gas heat capasity
%
% Output
% Ur        [=] J/m^2sK         Heat coefficient
% LAMBDAer  [=] J/msK           Effective radial conductivity


function [Ur,LAMBDAer]=heatcoef(Re,T,Y,VIS,CPgas)

global Dp RADIUSi RADIUSo LAMBDAst EPS RP LAMBDA Ncomp

P=1.0;
BETA=1.0;
LAMBDAs=0.243;
PHI=0.3;

%Calculates the gas heat conductivity
Tmatrix=[ones(RP,1) T T.^2 T.^3];	
for i=1:Ncomp
   LAMBDAcomp(:,i)=Tmatrix*LAMBDA(i,:)';				
end
LAMBDAg=diag(Y*LAMBDAcomp');

%Prandtl number
Pr = VIS.*CPgas./LAMBDAg;

%Radial effective static conduction
ALPHArv = (0.227/(1+EPS/(2*(1-EPS))*(1-P)/P)*(T/100).^3);
ALPHArs = 0.227*P/(2-P)*(T/100).^3;
LAMBDAer0 = LAMBDAg.*(EPS*(1 + BETA*Dp*ALPHArv./LAMBDAg) + ...
    BETA*(1-EPS)./(1./(1/PHI + ALPHArs*Dp./LAMBDAg) + 2/3*LAMBDAg/LAMBDAs));

%Effective radial conductivity
LAMBDAer = LAMBDAer0+0.14*LAMBDAg.*Re.*Pr;

%Heat transfer coefficient near the wall
ALPHAw0=8.694/(2*RADIUSi)^(4/3)*LAMBDAer0(RP);
ALPHAw=ALPHAw0+0.444*Re(RP)*Pr(RP)*LAMBDAg(RP)/Dp;

%Overall heat transfer coefficient
Ur=(RADIUSi*log(RADIUSo/RADIUSi)/LAMBDAst+1/ALPHAw)^(-1);
