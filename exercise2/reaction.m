% reaction
% This function calculates the reaction rates for all the components
% and the heat of the reaction in all the discretication points.
% Input:
% T         [=] K               Temperature
% Y         [=] -               Mol fraction
% P         [=] Pa              Total pressure
%
% Output
% Rcomp     [=] kmol/kg(cat)s   Reaction rate
% DELTAHr   [=] J/m^3s          Total reaction heat

function [Rcomp,DELTAHr] = reaction(T,Y,P)

global AJ ACTEN GASCONST AX ADENT RP RHOcat EPS MMASS ...
    ENT298 ENT948 Ncomp

%Reaction enthalpies 
for i=1:3
   DELTAH(:,i)=ENT298(i)+(T-298)/(948-298)*(ENT948(i)-ENT298(i));
end

%Partial pressures
for i=1:Ncomp
    Pcomp(:,i) = Y(:,i).*P/1e5;
end

DELTAHr = zeros(RP,1);
for i=1:RP
%Rate constant
  Krx  = AJ.*exp(-ACTEN./(GASCONST*T(i)))/3600;

%Adsorbtion constant
  Kads = AX.*exp(-ADENT./(GASCONST*T(i)));

%Equilibrium constants
  Keq(1) = 10^(-11650/T(i)+13.076);
  Keq(2) = 10^(1910/T(i)-1.784);
  Keq(3) = Keq(1)*Keq(2);

%Rx rates
  DENOM(i) = 1 + Kads(2)*Pcomp(i,2) + Kads(3)*Pcomp(i,4) ...
      + Kads(1)*Pcomp(i,1) + Kads(4)*Pcomp(i,5)/Pcomp(i,4);

  Rrx(1)= Krx(1)/(Pcomp(i,4))^2.5*(Pcomp(i,1)*Pcomp(i,5) ...
      - (Pcomp(i,4))^3*Pcomp(i,2)/Keq(1))/(DENOM(i))^2;
  Rrx(2)= Krx(2)/(Pcomp(i,4))*(Pcomp(i,2)*Pcomp(i,5) ...
      - (Pcomp(i,4))*Pcomp(i,3)/Keq(2))/(DENOM(i))^2;
  Rrx(3)= Krx(3)/(Pcomp(i,4))^3.5*(Pcomp(i,1)*(Pcomp(i,5))^2 ...
      - (Pcomp(i,4))^4*Pcomp(i,3)/Keq(3))/(DENOM(i))^2;

%Production/consumption rates of components
  Rcomp(i,1)=-(Rrx(1)+Rrx(3));
  Rcomp(i,2)=Rrx(1)-Rrx(2);
  Rcomp(i,3)=Rrx(2)+Rrx(3);
  Rcomp(i,4)=3.*Rrx(1)+Rrx(2)+4.*Rrx(3);
  Rcomp(i,5)=-(Rrx(1)+Rrx(2)+2.*Rrx(3));
  Rcomp(i,6)=0;

%Rx heat
  for j=1:3
    DELTAHr(i) = DELTAHr(i) + DELTAH(i,j)*Rrx(j)*RHOcat*(1-EPS);
  end
end
