function []=plottingFunc(start,r,z,y,ndisk)
var = y(:,start*ndisk+1:(start+1)*ndisk);
figure
mesh(r,z,var)
% axis([0 LENGTH 0 RADIUSi 0 max(T)])
grid on
xlabel('Radius [m]')
ylabel('Z [m]')
if(start==0)
    zlabel('CH4 mass fraction [-]');
end
if(start==1)
    zlabel('CO mass fraction [-]');
end
if(start==2)
    zlabel('CO2 mass fraction [-]');
end
if(start==3)
    zlabel('H2 mass fraction [-]');
end
if(start==4)
    zlabel('H2O mass fraction [-]');
end
if(start==5)
    zlabel('Temperature [K]');
end
if(start==6)
    zlabel('Uz [m/s]');
end

