function [CPgas]=cp(Y,T)
global Ncomp CP RP

%Calculates the gas heat conductivity
Tmatrix=[ones(RP,1) T T.^2 T.^3];	
for i=1:Ncomp
   CPcomp(:,i)=Tmatrix*CP(i,:)';				
end

for i=1:RP
    CPgas(i)= sum(CPcomp(i,:).* Y(i,:));
end

end