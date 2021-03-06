function [lambagas]=LAMBDAG(Y,T)
global Ncomp RP LAMBDA

%Calculates the gas heat conductivity
Tmatrix=[ones(RP,1) T T.^2 T.^3]; %[K]
for i=1:Ncomp
   LAMBDAcomp(:,i)=Tmatrix*LAMBDA(i,:)'; % [W/mK]		
end

for i=1:RP
    lambagas(i)= sum(LAMBDAcomp(i,:).* Y(i,:));
end

end