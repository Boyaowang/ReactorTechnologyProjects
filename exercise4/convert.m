% convert
% The function takes the component densities or the massfraction as
% input and converts it to mol fraction.
% Input:
% COMP      [=] -               Mass fraction or component density
%
% Output
% Y         [=] -               Mol fraction

function Y = convert(COMP)

global MMASS Ncomp

% Calculates the total density by summing up the component densities
RHOtot = 0;
for i=1:Ncomp
  RHOtot = RHOtot+COMP(i);
end

% Converts the input to molefraction by the following expression:
% Y=W(i)/MMASS(i)/SUM(W(i)/MMASS(i))
W   = 0;
SUM = 0;
for i=1:Ncomp
  W(i) = COMP(i)/RHOtot;
  SUM  = SUM + W(i)/MMASS(i);
end

for i=1:Ncomp
  Y(i) = W(i)/(MMASS(i)*SUM);
end
