% masscoef
% The function computes the molecular diffusivities and the mass transfer
% coefficients of the components.
% Input:
% Ptot      [=] Pa              Total pressure
% T         [=] K               Temperature
% viscg     [=] kg/ms           Gas viscosity
% rhog      [=] kg/m^3          Gas density
% vz        [=] m/s             Superficial velocity
%
% Output
% Dm        [=] m^2/s           Molecular diffusivity
% k         [=] m/s             Mass transfer coefficient

function [Dm,k]=masscoef(Ptot,T,rhog,vz,viscg)

global MMASS Ncomp mpart Dp sumny
for j=1:Ncomp
for i=1:Ncomp
% Calculation of the diffusion coefficients
    Ma(i)= 2*(1/MMASS(i)+1/MMASS(j))^(-1);
    Dm(:,i)=1e-4*0.00143*T.^1.75./...
        (Ptot*1e-5*sqrt(Ma(i))*(sumny(i)^(1/3)+sumny(j)^(1/3))^2);
    
% Calculation of the mass transfer coefficients    
    k(:,i) = vz.*1.17.*(Dp*vz.*rhog./viscg).^(-0.42)...
        .*(viscg./(rhog.*Dm(:,i))).^(-0.67);
end
end
