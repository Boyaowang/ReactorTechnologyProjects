%% heat capacity
function Cp = CpSO2(T)
Cp = 30.178+ T*42.452*10^(-3) - T*T*18.218*10^(-6);
end

function Cp = CpO2(T)
Cp = 23.995+ T*17.507*10^(-3) - T*T*6.628*10^(-6);
end

function Cp = CpSO3(T)
Cp = 35.634+ T*71.722*10^(-3) - T*T*31.539*10^(-6);
end

function Cp = CpN2(T)
Cp = 26.159+ T*6.615*10^(-3) - T*T*2.889^(-7);
end

%% Equilibrium constant
function Kp_ = Kp(T)
global R;
Kp_ = 3.142*10^(-3)*exp(98325/(R*T)-11.2); %[Pa^(-0.5)]
end

%% Equilibrium constant
function kRate = k(T)
kRate = 9.8692*10^(-3)*exp(-97782/T-110.1*log(T)+848.1); %[mol(SO2)/kg(cat)*Pa]
end

%% Reynolds number
function Rey = Rep(rhog,us,dp,mu)
Rey = rhog*us*dp/mu; %Calculating the friction factor;
end

%% Friction coefficient
function fri = f(epsilon, Reynolds)
fri = ((1-epsilon)/epsilon^3) * (1.75+ 4.2*Reynolds^(5/6)*(1-epsilon)/Reynolds);
end


