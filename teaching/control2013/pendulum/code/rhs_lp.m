% RHS of linear model
% no control, no noise
function dxdt = rhs_lp(t,x,M,m,l,g,k,c,I)

dxdt = zeros(4,1);

v1 = (M+m)/(I*(M+m)+l^2*m*M);
v2 = (I+l^2*m)/(I*(M+m)+l^2*m*M);

dxdt = [x(2);...
        (-k*v2*x(2) -m^2*l^2*g*v2*x(3)/(I + m*l^2) + m*l*c*v2*x(4)/(I + m*l^2));...
        x(4);...
        (m*l*k*v2*x(2)/(M+m) + m*l*g*v1*x(3) - c*v1*x(4))] ;

