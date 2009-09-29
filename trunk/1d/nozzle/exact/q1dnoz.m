clear all
close all

global gamma H0 frate L

gamma=1.4;

% Inflow conditions
m1=1.5; % mach number
r1=1.0; % density
p1=1.0;
c1=sqrt(gamma*p1/r1);
u1=m1*c1;
s1=p1/r1^gamma;

% Enthalpy: this is constant
H0 = c1^2/(gamma-1.0) + 0.5*u1^2;

% pressure ratio: Poutlet/Pinlet
prat=2.5;

% outflow pressure
p2=p1*prat;

% x location of nozzle inlet and outlet
x1 = 0.0; % inlet
a1 = nozarea(x1); % inlet area
x2 = L; % L is set when nozarea is called above

% flow rate: this is constant
frate = r1*u1*a1;

%x = linspace(0,4.5,30);
%for j=1:30
%r(j) = density(x(j), s1, r1);
%end
%
%load flow.dat;
%plot(x,r,'o',flow(:,1),flow(:,2))

% Find shock location xs
fun = @(x) shockfun(x, x2, r1, s1, p2);
xs = fminbnd(fun, x1, x2, optimset('TolX',1e-12))
fun(xs)

% Load numerical solution
load flow.dat

% plots
x = linspace(x1,x2,100);
[r u p m] = solution(xs, s1, r1, x);

figure(1)
plot(x,r,flow(:,1),flow(:,2),'o');
legend('Exact', 'Numerical')
title('Density')

figure(2)
plot(x,u,flow(:,1),flow(:,3),'o');
legend('Exact', 'Numerical')
title('Velocity')

figure(3)
plot(x,p,flow(:,1),flow(:,4),'o');
legend('Exact', 'Numerical')
title('Pressure')

figure(4)
plot(x,m,flow(:,1),flow(:,5),'o');
legend('Exact', 'Numerical')
title('Mach number')
