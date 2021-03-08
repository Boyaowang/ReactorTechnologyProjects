function dydr = calc(vars)
%%%%%%%%%%%%%%%%%%% Import constants:%%%%%%%%%%%%%%%%%%%%%%
global MM mpart GASCONST RHOcat rp EPS Ncomp MMASS Dp TEMPout rp pin

%%%%%%%%%%%% create variables vector for each z %%%%%%%%%%%%%%
r = (0:rp/mpart:rp);
wCH4 =vars(1:mpart);
wCO = vars(mpart+1:2*mpart) ;
wCO2 = vars(2*mpart+1:3*mpart);
wH2 = vars(3*mpart+1:4*mpart);
wH2O =vars(4*mpart+1:5*mpart);
wN2 = ones(mpart,1) - wCH4 - wCO -wH2 -wH2O;
T = vars(5*mpart+1:6*mpart);

%Calculation of rhog:
rhog = (pin*MM)./(GASCONST*T);

%%%%%%%%%%%%%%%%%%%%%%%% Parameter calculations %%%%%%%%%%%%%%%%%%%%
%Diffusivity
% Di=Dp*uin/(1.1 * 8 *(2- (1-Dp/rp)^2));

%Define mass and molar fractions:
Ymass =[wCH4 wCO wCO2 wH2 wH2O wN2];
Ymol = zeros(mpart, 6);
for i=1:mpart
    Ymol(i,:) =  convert(Ymass(i,:));
end

%Create the Matrix Rcomp and column vector for the enthalpy Note that these
%are GIven For the MOLAR VALUES, NOT MASS => need to convert to mass based
[Rcomp,DELTAHr] = reaction(T,Ymol,pin*ones(mpart,1));

%%%%%%%%%%%%%%%%%%%%%%%%% generate the solver %%%%%%%%%%%%%%%%%%%%
%Create 1st and 2nd derivative vector:
dwCH4dr = dss020(r(1),rp,mpart,wCH4,1)';
dwCH4dr2 = dss042(r(1),rp,mpart,wCH4,dwCH4dr,2,2)';

for i=2:mpart-1
    dwCH4dr2(i) = dwCH4dr2(i) +RHOcat.*Rcomp(i,1)./(rhog(i));
end

dwCOdr = dss020(r(1),rp,mpart,wCO,1)';
dwCOdr2 = dss042(r(1),rp,mpart,wCO,dwCOdr,2,2)';

for i=2:mpart-1
    dwCOdr2(i) = dwCOdr2(i) +RHOcat.*Rcomp(i,2)./(rhog(i));
end

dwCO2dr = dss020(r(1),rp,mpart,wCO2,1)';
dwCO2dr2 = dss042(r(1),rp,mpart,wCO2,dwCO2dr,2,2)';

for i=2:mpart-1
    dwCO2dr2(i) = dwCO2dr2(i) +RHOcat.*Rcomp(i,3)./(rhog(i));
end

dwH2dr = dss020(r(1),rp,mpart,wH2,1)';
dwH2dr2 = dss042(r(1),rp,mpart,wH2,dwH2dr,2,2)';

for i=2:mpart-1
    dwH2dr2(i) = dwH2dr2(i) +RHOcat.*Rcomp(i,4)./(rhog(i));
end

dwH2Odr = dss020(r(1),rp,mpart,wH2O,1)';
dwH2Odr2 = dss042(r(1),rp,mpart,wH2O,dwH2Odr,2,2)';

for i=2:mpart-1
    dwH2Odr2(i) = dwH2Odr2(i) +RHOcat.*Rcomp(i,5)./(rhog(i));
end

dTdr = dss020(r(1),rp,mpart,T,1)';
dTdr2 = dss042(r(1),rp,mpart,T,dTdr,2,2)';

for i=2:mpart-1
    dTdr2(i) = dTdr2(i) +r(i).^2.*DELTAHr(i);
end


% boundaries particle center 
dwCH4dr0 = dwCH4dr(1);
dwCOdr0 = dwCOdr(1);
dwCO2dr0 = dwCO2dr(1);
dwH2dr0 = dwH2dr(1);
dwH2Odr0 = dwH2Odr(1);
dTdr0 = dTdr(1);

% boundaries particle surface 
dwCH4drR = dwCH4dr(mpart);
dwCOdrR = dwCOdr(mpart);
dwCO2drR = dwCO2dr(mpart);
dwH2drR = dwH2dr(mpart);
dwH2OdrR = dwH2Odr(mpart);
dTdrR = dTdr(mpart);

dydr = [dwCH4dr0; dwCH4dr2(2:mpart-1); dwCH4drR;...
    dwCOdr0; dwCOdr2(2:mpart-1);dwCOdrR;...
    dwCO2dr0; dwCO2dr2(2:mpart-1); dwCO2drR;...
    dwH2dr0; dwH2dr2(2:mpart-1); dwH2drR;...
    dwH2Odr0; dwH2Odr2(2:mpart-1); dwH2OdrR;...
    dTdr0; dTdr2(2:mpart-1); dTdrR];
end
