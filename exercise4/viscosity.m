function VIS = viscosity(Y,T,np)
global S B Ncomp %viscosity coefficients

for i = 1:Ncomp
    muy(:,i) = (T.^1.5.*B(i))./(T+S(i));
end

%Viscosity:
for i=1:np
    VIS(i)= muy(i,:) * Y(i,:)';
end
end