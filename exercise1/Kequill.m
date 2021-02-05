
function K = Kp(T)
%Modified for R
R = 8314; %[J/kmole*K] should be global
K = 3.142*10^(-3)*exp(98325/(R*T)-11.2); %[Pa^(-0.5)]
end