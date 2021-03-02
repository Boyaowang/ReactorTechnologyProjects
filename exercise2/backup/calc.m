function dydz = calc(z,vari,par,r)

global MM RP GASCONST RHOcat RADIUSi EPS Ncomp MMASS Dp TEMPout

r0 = par(1); % Lower integration limit in r-direction
eta=par(2); % Efficency factor
n =RP; % number of descretization point
uin = par(3); % initial velocity profile

wCH4 =vari(1:n);
wCO = vari(n+1:2*n) ;
wCO2 = vari(2*n+1:3*n);
wH2 = vari(3*n+1:4*n);
wH2O =vari(4*n+1:5*n);
wN2 = ones(n,1) - wCH4 - wCO -wH2 -wH2O;
T = vari(5*n+1:6*n);
uz = vari(6*n+1:7*n);
ptot = vari(7*n+1);

Di=Dp*uin/(1.1 * 8 *(2- (1-Dp/RADIUSi)^2)); %to ask

Ymass =[wCH4 wCO wCO2 wH2 wH2O wN2];
Ymol = zeros(RP, 6);
for i=1:RP
    Ymol(i,:) =  convert(Ymass(i,:));
end
% Ymol = convert(Ymass);
%Gas heat capacity:
Cpg = cp(Ymol,T)'; %todo

%Create the Matrix Rcomp and column vector for the enthalpy Note that these
%are GIven For the MOLAR VALUES, NOT MASS
[Rcomp,DELTAHr] = reaction(T,Ymol,ptot*ones(n,1));

for i=1:Ncomp
    Rcomp(:,i) = Rcomp(:,i).*MMASS(i); %to ask!!!
end

%Gas conductivity
lambdar = LAMBDAG(Ymol,T)';

%Calculation of rhog:

rhog = (ptot*MM)./(GASCONST*T);

mu = viscosity(Ymol, T)';
Re =rhog.*uz.*Dp./mu; %todo

dwCH4dr = dss020(r0,RADIUSi,n,wCH4,1)';
dwCH4dr2 = dss042(r0,RADIUSi,n,wCH4,dwCH4dr,2,2)';

dwCOdr = dss020(r0,RADIUSi,n,wCO,1)';
dwCOdr2 = dss042(r0,RADIUSi,n,wCO,dwCOdr,2,2)';

dwCO2dr = dss020(r0,RADIUSi,n,wCO2,1)';
dwCO2dr2 = dss042(r0,RADIUSi,n,wCO2,dwCO2dr,2,2)';

dwH2dr = dss020(r0,RADIUSi,n,wH2,1)';
dwH2dr2 = dss042(r0,RADIUSi,n,wH2,dwH2dr,2,2)';

dwH2Odr = dss020(r0,RADIUSi,n,wH2O,1)';
dwH2Odr2 = dss042(r0,RADIUSi,n,wH2O,dwH2Odr,2,2)';

dTdr = dss020(r0,RADIUSi,n,T,1)';
dTdr2 = dss042(r0,RADIUSi,n,T,dTdr,2,2)';

duzdr =dss020(r0,RADIUSi,n,uz,-1)';
dpdz = ergun(rhog, uz, Re, r);

% transfer mole based to mass based Cpg and DELTAHr
Cpg = Cpg./MM;

[Ur,LAMBDAer]=heatcoef(Re,T,Ymol,mu,Cpg);

dTdz = (1./(rhog(2:n-1) .* uz(2:n-1) .* Cpg(2:n-1))).*...
       (LAMBDAer(2:n-1) .*(dTdr2(2:n-1) + 1./r(2:n-1) .* dTdr(2:n-1))...
       -DELTAHr(2:n-1).*eta);

resdTdz1 = dTdr(1);
resdTdz2 = dTdr(n)+Ur/LAMBDAer(RP)*(T(RP)-TEMPout); %% todo


% duzdz = -durdr(2:n-1) + ur(2:n-1) ...
%     .* (1./T(2:n-1) .* dTdr(2:n-1)- 1./r(2:n-1)) ...
%     + uz(2:n-1)./ptot(2:n-1) ...
%     .* (ptot(2:n-1)./T(2:n-1) .* dTdz - dpdz*ones(n-2,1));

duzdz = -uz(2:n-1)./ptot ...
    .* (dpdz*ones(n-2,1) - ptot./T(2:n-1) .* dTdz);

resduzdz1 =duzdr(1);
resduzdz2 =duzdr(n);

dwCH4dz = 1./uz(2:n-1) ...
    .* (Di*(1./r(2:n-1).*dwCH4dr(2:n-1) + dwCH4dr2(2:n-1))...
    - wCH4(2:n-1) .* duzdz ...
    + RHOcat*eta.*Rcomp(2:n-1,1)./(rhog(2:n-1)));
resdwCH4dr1 = dwCH4dr(1);
resdwCH4dr2 = dwCH4dr(n);

dwCOdz = 1./uz(2:n-1) ...
    .* (Di*(1./r(2:n-1).*dwCOdr(2:n-1) + dwCOdr2(2:n-1))...
    - wCO(2:n-1) .* duzdz ...
    + RHOcat*eta.*Rcomp(2:n-1,2)./(rhog(2:n-1)));
resdwCOdr1 = dwCOdr(1);
resdwCOdr2 = dwCOdr(n);

dwCO2dz = 1./uz(2:n-1) ...
    .* (Di*(1./r(2:n-1).*dwCO2dr(2:n-1) + dwCO2dr2(2:n-1))...
    - wCO2(2:n-1) .* duzdz...
    + RHOcat*eta.*Rcomp(2:n-1,3)./(rhog(2:n-1)));
resdwCO2dr1 = dwCO2dr(1);
resdwCO2dr2 = dwCO2dr(n);

dwH2dz = 1./uz(2:n-1) ...
    .* (Di*(1./r(2:n-1).*dwH2dr(2:n-1) + dwH2dr2(2:n-1))...
    - wH2(2:n-1) .* duzdz...
    + RHOcat*eta.*Rcomp(2:n-1,4)./(rhog(2:n-1)));
resdwH2dr1 = dwH2dr(1);
resdwH2dr2 = dwH2dr(n);

dwH2Odz = 1./uz(2:n-1) ...
    .* (Di*(1./r(2:n-1).*dwH2Odr(2:n-1) + dwH2Odr2(2:n-1))...
    - wH2O(2:n-1) .* duzdz...
    + RHOcat*eta.*Rcomp(2:n-1,5)./(rhog(2:n-1)));
resdwH2Odr1 = dwH2Odr(1);
resdwH2Odr2 = dwH2Odr(n);


dydz = [resdwCH4dr1;dwCH4dz;resdwCH4dr2;...
    resdwCOdr1;dwCOdz;resdwCOdr2;...
    resdwCO2dr1;dwCO2dz;resdwCO2dr2;...
    resdwH2dr1;dwH2dz;resdwH2dr2;...
    resdwH2Odr1;dwH2Odz;resdwH2Odr2;...
    resdTdz1;dTdz;resdTdz2;...
    resduzdz1;duzdz;resduzdz2;...
    dpdz];

end
