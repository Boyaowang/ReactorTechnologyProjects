function dydz = calcGas(z,vari,par,r)
%%%%%%%%%%%%%%%%%%% Import constants:%%%%%%%%%%%%%%%%%%%%%%
global RP GASCONST RHOcat RADIUSi EPS Ncomp MMASS Dp TEMPout mpart hp av
r0 = par(1); % Lower integration limit in r-direction
eta=par(2); % Efficency factor
uin = par(3); % initial velocity profile

%%%%%%%%%%%% Initialization for reactor %%%%%%%%%%%%%%
n = RP; % number of discretization point
wCH4 =vari(1:n)
wCO = vari(n+1:2*n) ;
wCO2 = vari(2*n+1:3*n);
wH2 = vari(3*n+1:4*n);
wH2O =vari(4*n+1:5*n);
wN2 = ones(n,1) - wCH4 - wCO -wH2 -wH2O-wCO2;
T = vari(5*n+1:6*n);
uz = vari(6*n+1:7*n);
ptot = vari(7*n+1);

for i=1:RP
    w_guess = [wCH4(i) wCO(i) wCO2(i) wH2(i) wH2O(i) wN2(i)];
    T_guess = T(i);
    uz_guess = uz(i);
    ptot_guess = ptot;
    yPellet(:,i) = solvePellet(w_guess, T_guess, ptot_guess, uz_guess);    
end

wCH4s = yPellet(mpart,:);
wCOs = yPellet(2*mpart,:) ;
wCO2s = yPellet(3*mpart,:);
wH2s = yPellet(4*mpart,:);
wH2Os =yPellet(5*mpart,:);
wN2s = 1 - wCH4s - wCO2s - wCOs -wH2s -wH2Os;

Ts = yPellet(6*mpart,:);

%%%%%%%%%%%%%%%%%%%%%%%% Parameter calculations %%%%%%%%%%%%%%%%%%%%
%Diffusivity
Di=Dp*uin/(1.1 * 8 *(2- (1-Dp/RADIUSi)^2));

%Define mass and molar fractions for gas phase:
Ymass =[wCH4 wCO wCO2 wH2 wH2O wN2];
Ymol = zeros(RP, 6);
for i=1:RP
    Ymol(i,:) =  convert(Ymass(i,:));
end

%Define mass and molar fractions for solid phase:
YmassS =[wCH4s' wCOs' wCO2s' wH2s' wH2Os' wN2s'];
YmolS = zeros(RP, 6);
for i=1:RP
    MM(i) = Ymass(i,:) * MMASS';
    YmolS(i,:) =  convert(YmassS(i,:));
end
MM = MM';
%Gas heat capacity:
Cpg = cp(Ymol,T)';
%Viscocity for gas:
mu = viscosity(Ymol, T, RP)';

%Viscocity for solid:
%VIS = viscosity(YmolS, Ts', mpart)';
for i=1:RP
    VIS(i) = viscosityS(YmolS(i,:),Ts(i));
end

%Calculation of rhog:
rhog = (ptot*MM)./(GASCONST*T);
%Reynold number :
Re =rhog.*uz.*Dp./mu;

% transfer mole based to mass based Cpg and DELTAHr
Cpg = Cpg./MM;
[Ur,LAMBDAer]=heatcoef(Re,T,Ymol,mu,Cpg);
%Mass transfer coef:
for i=1:RP
[Dim(i,:),k(i,:)]=masscoef(ptot,Ts(i),rhog(i),uz(i),VIS(i)...
    ,YmassS(i,:),MM(i));
end
%%%%%%%%%%%%%%%%%%%%%%%%% generate the solver %%%%%%%%%%%%%%%%%%%%
%Create 1st and 2nd derivative vector:
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
%%%%%%%%%%%%%%%%%%%%%%% %Pressure drop: %%%%%%%%%%%%%%%%%%%%%%%%%%
dpdz = ergun(rhog, uz, Re, r);

%%%%%%%%%%%%%%%%%%%%%%% %Temperature equation: %%%%%%%%%%%%%%%%%%%%%%%
dTdz = (1./(rhog(2:n-1) .* uz(2:n-1) .* Cpg(2:n-1).*EPS)).*...
       (LAMBDAer(2:n-1) .*(dTdr2(2:n-1) + 1./r(2:n-1) .* dTdr(2:n-1))...
       -hp*av*(T(2:n-1)-Ts(2:n-1)')); %%Ts,h,av
%BCs:
resdTdz1 = dTdr(1);
resdTdz2 = dTdr(n)+Ur/LAMBDAer(RP)*(T(RP)-TEMPout); %% todo

%%%%%%%%%%%%%%%%%%%%%%% %Velocity equation:%%%%%%%%%%%%%%%%%%%%%%%
%drho/dz
drhogdz = (MM(2:n-1)/GASCONST).*( (1./(T(2:n-1))./(dpdz*ones(n-2,1))) ...
                                - (ptot./(T(2:n-1).^2).*dTdz));

%drho/dr                            
drhogdr =  dss020(r0,RADIUSi,n,rhog,1)';

duzdz = - 1./rhog(2:n-1).*(drhogdz.*uz(2:n-1));
% ./EPS+...
%     +k((2:n-1),1) .* av./EPS.*(wCH4(2:n-1)-wCH4s(2:n-1)')...
%     +k((2:n-1),2) .* av./EPS.*(wCO(2:n-1)-wCOs(2:n-1)')...
%     +k((2:n-1),3) .* av./EPS.*(wCO2(2:n-1)-wCO2s(2:n-1)')...
%     +k((2:n-1),4) .* av./EPS.*(wH2(2:n-1)-wH2s(2:n-1)')...
%     +k((2:n-1),5) .* av./EPS.*(wH2O(2:n-1)-wH2Os(2:n-1)'));
% duzdz = -uz(2:n-1)./ptot ...
%     .* (dpdz*ones(n-2,1) - ptot./T(2:n-1) .* dTdz);
%BCs
resduzdz1 =duzdr(1);
resduzdz2 =duzdr(n);

%%%%%%%%%%%%%%%%%%%%%%% %CH4 - mass fraction:%%%%%%%%%%%%%%%%%%%%%%%
% dwCH4dz = 1./uz(2:n-1) ...
%     .* (Di/EPS*(1./r(2:n-1).*dwCH4dr(2:n-1) + dwCH4dr2(2:n-1))...
%     - wCH4(2:n-1) .* duzdz ...
%     -k((2:n-1),1) .* av./EPS.*(wCH4(2:n-1)-wCH4s(2:n-1)')); %% ki, wch4s
dwCH4dz = 1./(rhog(2:n-1).*uz(2:n-1))...
    .*(-rhog(2:n-1).*wCH4(2:n-1).*duzdz-wCH4(2:n-1).*uz(2:n-1).*drhogdz...
 +Di.*(drhogdr(2:n-1).*dwCH4dr(2:n-1)...
      +rhog(2:n-1)./r(2:n-1).*dwCH4dr(2:n-1)...
      +rhog(2:n-1).*dwCH4dr2(2:n-1))...
      -k((2:n-1),1) .* rhog(2:n-1).* av.*(wCH4(2:n-1)-wCH4s(2:n-1)'));
%BCs:
resdwCH4dr1 = dwCH4dr(1);
resdwCH4dr2 = dwCH4dr(n);

%%%%%%%%%%%%%%%%%%%%%%% %CO - mass fraction:%%%%%%%%%%%%%%%%%%%%%%%
% dwCOdz = 1./uz(2:n-1) ...
%     .* (Di/EPS*(1./r(2:n-1).*dwCOdr(2:n-1) + dwCOdr2(2:n-1))...
%     - wCO(2:n-1) .* duzdz ...
%     -k((2:n-1),2) .* av./EPS.*(wCO(2:n-1)-wCOs(2:n-1)'));

dwCOdz = 1./(rhog(2:n-1).*uz(2:n-1))...
    .*(-rhog(2:n-1).*wCO(2:n-1).*duzdz-wCO(2:n-1).*uz(2:n-1).*drhogdz...
 +Di.*(drhogdr(2:n-1).*dwCOdr(2:n-1)...
      +rhog(2:n-1)./r(2:n-1).*dwCOdr(2:n-1)...
      +rhog(2:n-1).*dwCOdr2(2:n-1))...
      -k((2:n-1),2) .* rhog(2:n-1).* av.*(wCO(2:n-1)-wCOs(2:n-1)'));
resdwCOdr1 = dwCOdr(1);
resdwCOdr2 = dwCOdr(n);

% %%%%%%%%%%%%%%%%%%%%%%% %CO2 - mass fraction:%%%%%%%%%%%%%%%%%%%%%%%
% dwCO2dz = 1./uz(2:n-1) ...
%     .* (Di/EPS*(1./r(2:n-1).*dwCO2dr(2:n-1) + dwCO2dr2(2:n-1))...
%     - wCO2(2:n-1) .* duzdz...
%     -k((2:n-1),3) .* av./EPS.*(wCO2(2:n-1)-wCO2s(2:n-1)'));
dwCO2dz = 1./(rhog(2:n-1).*uz(2:n-1))...
    .*(-rhog(2:n-1).*wCO2(2:n-1).*duzdz-wCO2(2:n-1).*uz(2:n-1).*drhogdz...
 +Di.*(drhogdr(2:n-1).*dwCO2dr(2:n-1)...
      +rhog(2:n-1)./r(2:n-1).*dwCO2dr(2:n-1)...
      +rhog(2:n-1).*dwCO2dr2(2:n-1))...
      -k((2:n-1),3) .* rhog(2:n-1).* av.*(wCO2(2:n-1)-wCO2s(2:n-1)'));

resdwCO2dr1 = dwCO2dr(1);
resdwCO2dr2 = dwCO2dr(n);

%%%%%%%%%%%%%%%%%%%%%%%% H2 - mass fraction:%%%%%%%%%%%%%%%%%%%%%%%
% dwH2dz = 1./uz(2:n-1) ...
%     .* (Di/EPS*(1./r(2:n-1).*dwH2dr(2:n-1) + dwH2dr2(2:n-1))...
%     - wH2(2:n-1) .* duzdz...
%     -k((2:n-1),4) .* av./EPS.*(wH2(2:n-1)-wH2s(2:n-1)'));
dwH2dz = 1./(rhog(2:n-1).*uz(2:n-1))...
    .*(-rhog(2:n-1).*wH2(2:n-1).*duzdz-wH2(2:n-1).*uz(2:n-1).*drhogdz...
 +Di.*(drhogdr(2:n-1).*dwH2dr(2:n-1)...
      +rhog(2:n-1)./r(2:n-1).*dwH2dr(2:n-1)...
      +rhog(2:n-1).*dwH2dr2(2:n-1))...
      -k((2:n-1),4) .* rhog(2:n-1).* av.*(wH2(2:n-1)-wH2s(2:n-1)'));
resdwH2dr1 = dwH2dr(1);
resdwH2dr2 = dwH2dr(n);

%%%%%%%%%%%%%%%%%%%%%%%% H2O - mass fraction:%%%%%%%%%%%%%%%%%%%%%%%
% dwH2Odz = 1./uz(2:n-1) ...
%     .* (Di/EPS*(1./r(2:n-1).*dwH2Odr(2:n-1) + dwH2Odr2(2:n-1))...
%     - wH2O(2:n-1) .* duzdz...
%     -k((2:n-1),5) .* av./EPS.*(wH2O(2:n-1)-wH2Os(2:n-1)'));
dwH2Odz = 1./(rhog(2:n-1).*uz(2:n-1))...
    .*(-rhog(2:n-1).*wH2O(2:n-1).*duzdz-wH2O(2:n-1).*uz(2:n-1).*drhogdz...
 +Di.*(drhogdr(2:n-1).*dwH2Odr(2:n-1)...
      +rhog(2:n-1)./r(2:n-1).*dwH2Odr(2:n-1)...
      +rhog(2:n-1).*dwH2Odr2(2:n-1))...
      -k((2:n-1),5) .* rhog(2:n-1).* av.*(wH2O(2:n-1)-wH2Os(2:n-1)'));

resdwH2Odr1 = dwH2Odr(1);
resdwH2Odr2 = dwH2Odr(n);
% (dwCH4dz + dwCOdz + dwCO2dz + dwH2Odz + dwH2dz);
%%%%%%%%%%%%%%%%%%%%%%%%Define the equations system:%%%%%%%%%%%%%%%%%%%%%%
dydz = [resdwCH4dr1;dwCH4dz;resdwCH4dr2;...
    resdwCOdr1;dwCOdz;resdwCOdr2;...
    resdwCO2dr1;dwCO2dz;resdwCO2dr2;...
    resdwH2dr1;dwH2dz;resdwH2dr2;...
    resdwH2Odr1;dwH2Odz;resdwH2Odr2;...
    resdTdz1;dTdz;resdTdz2;...
    resduzdz1;duzdz;resduzdz2;...
    dpdz];

end
