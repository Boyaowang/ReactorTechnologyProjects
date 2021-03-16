clc;
close all;
clear all;
%%%%%%%%%%%%%%%%%%% Loading the constants %%%%%%%%%%%%%%%%
constant
eta = 0.001; % efficiency factor

%%%%%%%%%%%%%%%%% Sizing the r- and z-direction %%%%%%%%%%%%%
% r-direction:
r0 = 0; % Lower integration limit in r-direction
rn = rp; % Upper integration limit in r-direction
r = (r0:(rn-r0)/(mpart-1):rn)';

%%%%%%%%%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%
wCH40 = FRACin(1) .*ones(mpart,1);
wCO0 = FRACin(2) .*ones(mpart,1);
wCO20 = FRACin(3) .*ones(mpart,1);
wH20 = FRACin(4) .*ones(mpart,1);
wH2O0 = FRACin(5) .*ones(mpart,1);
wN2 = FRACin(6) .*ones(mpart,1);
T0 = Tin*ones(mpart,1);

%initial guess vectors:
init = [wCH40; wCO0; wCO20; wH20; wH2O0; T0];
%parameters vectors:
par = [r0 eta uin];
%%%%%%%%%%%%%%%%%%%%% non-linear system solver %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create ODE solver:
options = optimoptions('fsolve','Display','none','PlotFcn',@optimplotfirstorderopt);

func = @calc;
y = fsolve(@calc, init, options);



%Define the result vector:
n = mpart;
wCH4 = y(1:mpart);
wCO = y(mpart+1:2*mpart) ;
wCO2 = y(2*mpart+1:3*mpart);
wH2 = y(3*mpart+1:4*mpart);
wH2O =y(4*mpart+1:5*mpart);
wN2 = 1 - wCH4 - wCO2 - wCO -wH2 -wH2O;
T = y(5*mpart+1:6*mpart);

%%%%%%%%%%%%%%%%%% Ploting the results (species fraction)%%%%%%%%%%%%%%%%%%
figure()
for i=1:Ncomp-1
    m = 3;
    n = 2;
    subplot(m,n,i);
    plot(r, y((i-1)*mpart+1:i*mpart));

    xlabel('R [m]') 
    ylabel('y [-]')
    
if(i==1)
    title('CH4 mass fraction [-]');
end
if(i==2)
    title('CO mass fraction [-]');
end
if(i==3)
    title('CO2 mass fraction [-]');
end
if(i==4)
    title('H2 mass fraction [-]');
end
if(i==5)
    title('H2O mass fraction [-]');
end

end
subplot(m,n,6);
plot(r, wN2);
xlabel('R [m]') 
ylabel('y [-]')
title('N2 mass fraction [-]');
figure()
plot(r, T)
xlabel('r [m]') 
ylabel('T [K]')
title('Temperature [K]');

figure()
plot(r, wCH4);
xlabel('R [m]')
ylabel('y [-]')
title(['CH4 mass fraction (rp = ',num2str(rp),') [-]']);
