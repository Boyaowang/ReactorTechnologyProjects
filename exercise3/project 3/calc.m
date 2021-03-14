function dydr = calc(vars)
%%%%%%%%%%%%%%%%%%% Import constants:%%%%%%%%%%%%%%%%%%%%%%
global MM mpart GASCONST RHOcat rp EPS Ncomp MMASS Dp TEMPout rp pin uin...
    Tin FRACin

%%%%%%%%%%%% create variables vector for each z %%%%%%%%%%%%%%
r = linspace(0,rp,mpart);
wCH4 =vars(1:mpart);
wCO = vars(mpart+1:2*mpart) ;
wCO2 = vars(2*mpart+1:3*mpart);
wH2 = vars(3*mpart+1:4*mpart);
wH2O =vars(4*mpart+1:5*mpart);
wN2 = ones(mpart,1) - wCH4 - wCO -wH2 -wH2O - wCO2;
T = vars(5*mpart+1:6*mpart);

%Calculation of rhog:
rhog = (pin*MM)./(GASCONST*T);

%%%%%%%%%%%%%%%%%%%%%%%% Parameter calculations %%%%%%%%%%%%%%%%%%%%

%Define mass and molar fractions:
Ymass =[wCH4 wCO wCO2 wH2 wH2O wN2];
Ymol = zeros(mpart, 6);
for i=1:mpart
    Ymol(i,:) =  convert(Ymass(i,:));
end

%Viscocity:
VIS = viscosity(Ymol,T)';

%Diffusivity
for i = 1:mpart
    [Dim(i,:),k(i,:)]=masscoef(pin,T(i),rhog(i),uin,VIS(mpart),Ymass(i,:));
end

%Create the Matrix Rcomp and column vector for the enthalpy Note that these
%are GIven For the MOLAR VALUES, NOT MASS => need to convert to mass based
[Rcomp,DELTAHr] = reaction(T,Ymol,pin*ones(mpart,1));
for i=1:Ncomp % add (boyao)
    Rcomp(:,i) = Rcomp(:,i).*MMASS(i); %convert to mass based
end
%%%%%%%%%%%%%%%%%%%%%%%%% generate the solver %%%%%%%%%%%%%%%%%%%%

%Create 1st and 2nd derivative vector  +  residue vectors:

%Mass equation:

%%%%%%%% CH4 %%%%%%%%%%%%%%
dwCH4dr = dss020(r(1),rp,mpart,wCH4,-1)';
dwCH4dr2 = dss042(r(1),rp,mpart,wCH4,dwCH4dr,2,2)';
F_CH4 = rhog.*Dim(:,1).*(2.*dwCH4dr+r'.*dwCH4dr2) + RHOcat.*Rcomp(:,1).*(1-EPS);
% d1dr_CH4 = dss020(r(1),rp,mpart,wCH4,1)';
% dwCH4dr = -((r.^2)'.*(rhog.*Dim(:,1))).*dss020(r(1),rp,mpart,wCH4,1)';
% dwCH4dr2 = dss042(r(1),rp,mpart,wCH4,dwCH4dr,2,2)';
% F_CH4 = (-1./r').*dwCH4dr2 +RHOcat.*Rcomp(:,1).*(1-EPS); % add (boyao)

%%%%%%%% CO %%%%%%%%%%%%%%
dwCOdr = dss020(r(1),rp,mpart,wCO,-1)';
dwCOdr2 = dss042(r(1),rp,mpart,wCO,dwCOdr,2,2)';
F_CO = rhog.*Dim(:,2).*(2*dwCOdr+r'.*dwCOdr2) + RHOcat.*Rcomp(:,2).*(1-EPS);
% d1dr_CO = dss020(r(1),rp,mpart,wCO,1)';
% dwCOdr = -((r.^2)'.*(rhog.*Dim(:,2))).*dss020(r(1),rp,mpart,wCO,1)';
% dwCOdr2 = dss042(r(1),rp,mpart,wCO,dwCOdr,2,2)';
% F_CO = (-1./r').*dwCOdr2 +RHOcat.*Rcomp(:,2).*(1-EPS);

%%%%%%%% CO2 %%%%%%%%%%%%%%
% d1dr_CO2 = dss020(r(1),rp,mpart,wCO2,1)';
% dwCO2dr = -((r.^2)'.*(rhog.*Dim(:,3))).*dss020(r(1),rp,mpart,wCO2,1)';
% dwCO2dr2 = dss042(r(1),rp,mpart,wCO2,dwCO2dr,2,2)';
% F_CO2 = (-1./r').*dwCO2dr2 +RHOcat.*Rcomp(:,3).*(1-EPS);
dwCO2dr = dss020(r(1),rp,mpart,wCO2,-1)';
dwCO2dr2 = dss042(r(1),rp,mpart,wCO2,dwCO2dr,2,2)';
F_CO2 = rhog.*Dim(:,3).*(2*dwCO2dr+r'.*dwCO2dr2) + RHOcat.*Rcomp(:,3).*(1-EPS);

%%%%%%%% H2 %%%%%%%%%%%%%%
dwH2dr = dss020(r(1),rp,mpart,wH2,-1)';
dwH2dr2 = dss042(r(1),rp,mpart,wH2,dwH2dr,2,2)';
F_H2 = rhog.*Dim(:,4).*(2*dwH2dr+r'.*dwH2dr2) + RHOcat.*Rcomp(:,4).*(1-EPS);

%%%%%%% H2O %%%%%%%%%%%%%%
dwH2Odr = dss020(r(1),rp,mpart,wH2O,-1)';
dwH2Odr2 = dss042(r(1),rp,mpart,wH2O,dwH2Odr,2,2)';
F_H2O = rhog.*Dim(:,5).*(2.*dwH2Odr+r'.*dwH2Odr2) + RHOcat.*Rcomp(:,5).*(1-EPS);


%Temperature equation:
LAMBDA = 50;                                                    %Unfinished
% d1dr_T = dss020(r(1),rp,mpart,T,1)';
% dTdr = (r.^2)'.*(-LAMBDA*dss020(r(1),rp,mpart,T,1)');
% dTdr2 = (-1./r.^2)'.*dss042(r(1),rp,mpart,T,dTdr,2,2)';
% F_T = dTdr2 + RHOcat.*DELTAHr;
dTdr = dss020(r(1),rp,mpart,T,-1)';
dTdr2 = dss042(r(1),rp,mpart,T,dTdr,2,2)';
F_T = LAMBDA.*(2*dTdr + r'.*dTdr2) + r'.*DELTAHr;


% boundaries particle center 
F_CH4(1) = dwCH4dr(1);
F_CO(1) = dwCOdr(1);
F_CO2(1) = dwCO2dr(1);
F_H2(1) = dwH2dr(1);
F_H2O(1) = dwH2Odr(1);
F_T(1) = dTdr(1);
%%%%%%%%%%%%%%%%%%%%%%% Unfinished %%%%%%%%%%%%%%

k_CH4 = k(mpart,1);
wbulk_CH4 =FRACin(1);
k_CO = k(mpart,2);
wbulk_CO = FRACin(2);
k_CO2 = k(mpart,3);
wbulk_CO2 = FRACin(3);
k_H2 = k(mpart,4);
wbulk_H2 = FRACin(4);
k_H2O = k(mpart,5);
wbulk_H2O = FRACin(5);
%Heat transfer coefficient
h = 30000; %[W/m2K]
%Particle conductivity 
LAMBDA_b = 50; %[W/mK]
Tbulk = Tin;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% boundaries particle surface 
F_CH4(mpart) = dwCH4dr(mpart) + (k_CH4/Dim(mpart,1))*(wCH4(mpart) - wbulk_CH4);
F_CO(mpart) = dwCOdr(mpart) + (k_CO/Dim(mpart,2))*(wCO(mpart) - wbulk_CO);
F_CO2(mpart) = dwCO2dr(mpart) + (k_CO2/Dim(mpart,3))*(wCO2(mpart) - wbulk_CO2);
F_H2(mpart) = dwH2dr(mpart) + (k_H2/Dim(mpart,4))*(wH2(mpart) - wbulk_H2);
F_H2O(mpart) = dwH2Odr(mpart) + (k_H2O/Dim(mpart,5))*(wH2O(mpart) - wbulk_H2O);
F_T(mpart) = dTdr(mpart) + h/LAMBDA_b*(T(mpart) - Tbulk);

dydr = [F_CH4;F_CO;F_CO2;F_H2;F_H2O;F_T];
end
