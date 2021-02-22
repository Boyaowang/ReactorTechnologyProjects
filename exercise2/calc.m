function dydz = calc(z,vari,r)

Rt = 0.102; % tube radius [m]
r0 = 0; % Lower integration limit in r-direction
rn = Rt; % Upper integration limit in r-direction
ndisk = 50; % Number of discretization points
r = (r0:(rn-r0)/(ndisk-1):rn)';
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
T = vari(n+1:2*n) ;
us = vari(2*n+1:3*n);
ptot = vari(3*n+1:4*n);

rhog = 1* ones(n,1);


dwCH4dr = dss020(zstart,zend,n,wCH4,-1)';
dwCH4dr2 = dss042(zstart,zend,n,wCH4,dwCH4dr,2,2)';
dTdr = dss020(zstart,zend,n,T,-1)';
dTdr2 = dss042(zstart,zend,n,T,dTdr,2,2)';

dpdz = ergun(rhog, us, Re, r);
dTdz = (1./(rhog(2:n-1) .* us(2:n-1) .* Cpg)).*...
       (lambdar *(dTdr2(2:n-1) + 1/r((2:n-1)) .* dTdr((2:n-1))));
resdTdz1 = dTdr(1);
resdTdz2 = dTdr(n);

dusdz = -us./ptot * (dpdz - ptot./T .* dTdz);
dwCH4dz = 1/us * (Di/EPS*(1/r(2:n-1)*dwCH4dr(2:n-1) + dwCH4dr2(2:n-1))...
                  - wCH4(2:n-1) * dusdz(2:n-1));
resdwCH4dz1 = dwCH4dr(1);
resdwCH4dz2 = dwCH4dr(n)-100;

dydz = [resdwCH4dz1;dwCH4dz;resdwCH4dz2;resdTdz1;dTdz;resdTdz2;dpdz;dusdz];

end
