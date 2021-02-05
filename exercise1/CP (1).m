
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


