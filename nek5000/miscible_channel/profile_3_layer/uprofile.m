% Compute velocity profile for 3-layer miscible channel flow
% concentration profile is in c0.m
m = 10;
P = -1.0;

x0=0.5;
x1=1.0;

c = chebfun(@c0, [x0,x1], 'splitting', 'on');
mu = exp(c*log(m));
dmu= diff(mu);

L = chebop(x0,x1);
L.op = @(x,u) mu.*diff(u,2)+ dmu.*diff(u);
L.lbc = @(u) diff(u);
L.rbc = @(u) u;
u = L \ P;
plot(u)

% calculate flow rate
Q = 2*sum(u,x0,x1)

% normalize so that average velocity is one
u = u/Q;
plot(u)
xlabel('y')
ylabel('x velocity')

% Again calculate to check flow rate is one
Q = 2*sum(u,x0,x1);
assert(abs(Q-1.0) < 1.0e-12)

% To compute velocity at any point y
%   u(y)

nelem = 24
nx1   = 10
y = grid(nelem, nx1); % This y is in [0,0.5]
y = 1 - flipud(y);
uxdata = u(y);
cdata  = c0(y);

uxdata = [flipud(uxdata); uxdata];
cdata  = [flipud(cdata) ; cdata];
y      = [ 1-flipud(y)  ; y];

data = [y, uxdata, cdata];
save('3layer.dat','data','-ascii','-double')
