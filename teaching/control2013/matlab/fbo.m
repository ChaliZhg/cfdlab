% Nonlinear model
% no control, no noise
function dxdt = fbo(t,x,M,m,l,g,k,d,I)

dxdt = zeros(4,1);

dxdt = [ x(2);...
         ((l*m*d)/(I+l^2*m)*cos(x(3))*x(4) + l*m*x(4)^2*sin(x(3)) - (l*m)^2/(I+l^2*m)*g*cos(x(3))*sin(x(3)) - k*x(2) )/( (M+m) - (l*m)^2/(I+l^2*m)*cos(x(3))^2 );...
         x(4);...
        ( l*m*g*sin(x(3)) - (l*m)^2/(M+m)*x(4)^2*sin(x(3))*cos(x(3)) -d*x(4)  + k*l*m/(M+m)*cos(x(3))*x(2) )/( I + l^2*m - (l*m)^2/(M+m)*cos(x(3)^3)) ] ;
