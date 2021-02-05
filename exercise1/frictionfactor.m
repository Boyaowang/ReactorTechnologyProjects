
function Rey = Rep(rhog,us,dp,mu)
Rey = rhog*us*dp/mu; %Calculating the friction factor;
end
function fri = f(epsilon, Reynolds)
fri = ((1-epsilon)/epsilon^3) * (1.75+ 4.2*Reynolds^(5/6)*(1-epsilon)/Reynolds);
end
