function VIS = viscosityS(Y,T)
global S B Ncomp %viscosity coefficients

for i = 1:Ncomp
    muy(i) = (T^1.5*B(i))/(T+S(i));
end

VIS = 0;
for i = 1:Ncomp
    VIS = VIS+muy(i)*Y(i);
end

end