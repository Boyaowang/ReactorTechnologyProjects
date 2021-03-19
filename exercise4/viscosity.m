function VIS = viscosity(Y,T)
global S B Ncomp mpart %viscosity coefficients

for i = 1:Ncomp
    muy(:,i) = (T.^1.5.*B(i))./(T+S(i));
end

%Viscosity:
for i=1:mpart
    VIS(i)= muy(i,:) * Y(i,:)';
end
end