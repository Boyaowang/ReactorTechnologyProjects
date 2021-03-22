clc
fun = @root2d;
x0 = [0,0];
a = 5;
x = fsolve(@(x0)fun(x0,a),x0)

function F = root2d(x,a)
F(1) = exp(-exp(-(x(1)+x(2)))) - x(2)*(1+x(1)^2);
F(2) = x(1)*cos(x(2)) + x(2)*sin(x(1)) - 0.5;
a
end