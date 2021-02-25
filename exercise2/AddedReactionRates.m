function dydz = calc(z,vari,par,r)

Rt = 0.102; % tube radius [m]
r0 = 0; % Lower integration limit in r-direction
rn = Rt; % Upper integration limit in r-direction
ndisk = 50; % Number of discretization points
% r = (r0:(rn-r0)/(ndisk-1):rn)';
n =50;
zstart = 0; %[m]
zend = 6.096; %[m]
Cpg = 1000;
lambdar = -50;
EPS=0.5;
Di=0.01;
rhob=2000 ;
eta=0.5;
Re = 1000* ones(n,1);

wCH4 =vari(1:n);
wCO = vari(n+1:2*n) ;
wCO2 = vari(2*n+1:3*n);
wH2 = vari(3*n+1:4*n);
wH2O =vari(4*n+1:5*n);
T = vari(5*n+1:6*n) ;
uz = vari(6*n+1:7*n);
ptot = vari(7*n+1);

rhog = 1* ones(n,1);

%CHANGE: reating a matrix of the molar fractions
COMP = [wCH4 wCO wC02 wH2 wH20 wN2];
Y = convert(COMP);
%Create the Matrix Rcomp and column vector for the enthalpy Note that these
%are GIven For the MOLAR VALUES, NOT MASS
[Rcomp,DELTAHr] = reactions(T,Y,ptot);

dwCH4dr = dss020(r0,Rt,n,wCH4,-1)';
dwCH4dr2 = dss042(r0,Rt,n,wCH4,dwCH4dr,2,2)';

dwCOdr = dss020(r0,Rt,n,wCO,-1)';
dwCOdr2 = dss042(r0,Rt,n,wCO,dwCOdr,2,2)';

dwCO2dr = dss020(r0,Rt,n,wCO2,-1)';
dwCO2dr2 = dss042(r0,Rt,n,wCO2,dwCO2dr,2,2)';

dwH2dr = dss020(r0,Rt,n,wH2,-1)';
dwH2dr2 = dss042(r0,Rt,n,wH2,dwH2dr,2,2)';

dwH2Odr = dss020(r0,Rt,n,wH2O,-1)';
dwH2Odr2 = dss042(r0,Rt,n,wH2O,dwH2Odr,2,2)';

dTdr = dss020(r0,Rt,n,T,-1)';
dTdr2 = dss042(r0,Rt,n,T,dTdr,2,2)';

duzdr =dss020(r0,Rt,n,uz,-1)';
dpdz = ergun(rhog, uz, Re, r);

%CHANGE: Summed the reaction terms using the function sum
dTdz = (1./(rhog(2:n-1) .* uz(2:n-1) .* Cpg)).*...
       (lambdar .*(dTdr2(2:n-1) + 1./r((2:n-1)) .* dTdr((2:n-1))))...
       + (-DELTAHr) *0.001;
   %%todo
resdTdz1 = dTdr(1);
resdTdz2 = dTdr(n); %% todo


% duzdz = -durdr(2:n-1) + ur(2:n-1) ...
%     .* (1./T(2:n-1) .* dTdr(2:n-1)- 1./r(2:n-1)) ...
%     + uz(2:n-1)./ptot(2:n-1) ...
%     .* (ptot(2:n-1)./T(2:n-1) .* dTdz - dpdz*ones(n-2,1));

duzdz = -uz(2:n-1)./ptot ...
    .* (dpdz*ones(n-2,1) - ptot./T(2:n-1) .* dTdz);

resduzdz1 =duzdr(1);
resduzdz2 =duzdr(n);

dwCH4dz = 1./uz(2:n-1) ...
    .* (Di/EPS*(1./r(2:n-1).*dwCH4dr(2:n-1) + dwCH4dr2(2:n-1))...
    - wCH4(2:n-1) .* duzdz) + RHOcat*0.001*(1-EPS).*Rcomp(:,1);
resdwCH4dr1 = dwCH4dr(1);
resdwCH4dr2 = dwCH4dr(n);

dwCOdz = 1./uz(2:n-1) ...
    .* (Di/EPS*(1./r(2:n-1).*dwCOdr(2:n-1) + dwCOdr2(2:n-1))...
    - wCO(2:n-1) .* duzdz)+ RHOcat*0.001*(1-EPS).*Rcomp(:,2);
resdwCOdr1 = dwCOdr(1);
resdwCOdr2 = dwCOdr(n);

dwCO2dz = 1./uz(2:n-1) ...
    .* (Di/EPS*(1./r(2:n-1).*dwCO2dr(2:n-1) + dwCO2dr2(2:n-1))...
    - wCO2(2:n-1) .* duzdz) + RHOcat*0.001*(1-EPS).*Rcomp(:,3);
resdwCO2dr1 = dwCO2dr(1);
resdwCO2dr2 = dwCO2dr(n);

dwH2dz = 1./uz(2:n-1) ...
    .* (Di/EPS*(1./r(2:n-1).*dwH2dr(2:n-1) + dwH2dr2(2:n-1))...
    - wH2(2:n-1) .* duzdz)+ RHOcat*0.001*(1-EPS).*Rcomp(:,4);
resdwH2dr1 = dwH2dr(1);
resdwH2dr2 = dwH2dr(n);

dwH2Odz = 1./uz(2:n-1) ...
    .* (Di/EPS*(1./r(2:n-1).*dwH2Odr(2:n-1) + dwH2Odr2(2:n-1))...
    - wH2O(2:n-1) .* duzdz) + RHOcat*0.001*(1-EPS).*Rcomp(:,5);
resdwH2Odr1 = dwH2Odr(1);
resdwH2Odr2 = dwH2Odr(n);

%What about N2?


dydz = [resdwCH4dr1;dwCH4dz;resdwCH4dr2;...
    resdwCOdr1;dwCOdz;resdwCOdr2;...
    resdwCO2dr1;dwCO2dz;resdwCO2dr2;...
    resdwH2dr1;dwH2dz;resdwH2dr2;...
    resdwH2Odr1;dwH2Odz;resdwH2Odr2;...
    resdTdz1;dTdz;resdTdz2;...
    resduzdz1;duzdz;resduzdz2;...
    dpdz];

end