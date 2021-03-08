function dydr = calc(vars)
%%%%%%%%%%%%%%%%%%% Import constants:%%%%%%%%%%%%%%%%%%%%%%
global MM RP GASCONST RHOcat RADIUSi EPS Ncomp MMASS Dp TEMPout rp

%%%%%%%%%%%% create variables vector for each z %%%%%%%%%%%%%%
r = (0:rp/RP:rp);
wCH4 =vars(1:RP);
wCO = vars(RP+1:2*RP) ;
wCO2 = vars(2*RP+1:3*RP);
wH2 = vars(3*RP+1:4*RP);
wH2O =vars(4*RP+1:5*RP);
wN2 = ones(RP,1) - wCH4 - wCO -wH2 -wH2O;
T = vars(5*RP+1:6*RP);

%%%%%%%%%%%%%%%%%%%%%%%% Parameter calculations %%%%%%%%%%%%%%%%%%%%
%Diffusivity
% Di=Dp*uin/(1.1 * 8 *(2- (1-Dp/RADIUSi)^2));

%Define mass and molar fractions:
Ymass =[wCH4 wCO wCO2 wH2 wH2O wN2];
Ymol = zeros(RP, 6);
for i=1:RP
    Ymol(i,:) =  convert(Ymass(i,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%% generate the solver %%%%%%%%%%%%%%%%%%%%
%Create 1st and 2nd derivative vector:
dwCH4dr = dss020(r(1),RADIUSi,RP,wCH4,1)';
dwCH4dr2 = dss042(r(1),RADIUSi,RP,wCH4,dwCH4dr,2,2)';

dwCOdr = dss020(r(1),RADIUSi,RP,wCO,1)';
dwCOdr2 = dss042(r(1),RADIUSi,RP,wCO,dwCOdr,2,2)';

dwCO2dr = dss020(r(1),RADIUSi,RP,wCO2,1)';
dwCO2dr2 = dss042(r(1),RADIUSi,RP,wCO2,dwCO2dr,2,2)';

dwH2dr = dss020(r(1),RADIUSi,RP,wH2,1)';
dwH2dr2 = dss042(r(1),RADIUSi,RP,wH2,dwH2dr,2,2)';

dwH2Odr = dss020(r(1),RADIUSi,RP,wH2O,1)';
dwH2Odr2 = dss042(r(1),RADIUSi,RP,wH2O,dwH2Odr,2,2)';

dTdr = dss020(r(1),RADIUSi,RP,T,1)';
dTdr2 = dss042(r(1),RADIUSi,RP,T,dTdr,2,2)';

% boundaries particle center 
dwCH4dr0 = dwCH4dr(1);
dwCOdr0 = dwCOdr(1);
dwCO2dr0 = dwCO2dr(1);
dwH2dr0 = dwH2dr(1);
dwH2Odr0 = dwH2Odr(1);
dTdr0 = dTdr(1);

% boundaries particle surface 
dwCH4drR = dwCH4dr(RP);
dwCOdrR = dwCOdr(RP);
dwCO2drR = dwCO2dr(RP);
dwH2drR = dwH2dr(RP);
dwH2OdrR = dwH2Odr(RP);
dTdrR = dTdr(RP);

dydr = [dwCH4dr0; dwCH4dr2(2:RP-1); dwCH4drR;...
    dwCOdr0; dwCOdr2(2:RP-1);dwCOdrR;...
    dwCO2dr0; dwCO2dr2(2:RP-1); dwCO2drR;...
    dwH2dr0; dwH2dr2(2:RP-1); dwH2drR;...
    dwH2Odr0; dwH2Odr2(2:RP-1); dwH2OdrR;...
    dTdr0; dTdr2(2:RP-1); dTdrR];
end
