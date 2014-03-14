% Construct chebyshev interpolation and differentiation
% Run xconst.m before running this program
% Creates baseflow.mat file
clear all

ny1 = 10;

load xconst.mat

yd    = data2(:, 2);
ud    = data2(:, 3);
uyd   = data2(:, 4);
uyyd  = data2(:, 5);
uyyyd = data2(:, 6);
Td    = data2(:, 7);
Tyd   = data2(:, 8);
Tyyd  = data2(:, 9);
cd    = data2(:,10);
cyd   = data2(:,11);
cyyd  = data2(:,12);

u    = @(x) intgll(ny1, yd, ud,    x);
du   = @(x) intgll(ny1, yd, uyd,   x);
ddu  = @(x) intgll(ny1, yd, uyyd,  x);
dddu = @(x) intgll(ny1, yd, uyyyd, x);
T    = @(x) intgll(ny1, yd, Td,    x);
dT   = @(x) intgll(ny1, yd, Tyd,   x);
ddT  = @(x) intgll(ny1, yd, Tyyd,  x);
c    = @(x) intgll(ny1, yd, cd,    x);
dc   = @(x) intgll(ny1, yd, cyd,   x);
ddc  = @(x) intgll(ny1, yd, cyyd,  x);

%x = linspace(0,1,1000);
%plot(x,u(x))

uf = chebfun(u, [0,1], 'splitting', 'on');
duf = chebfun(du, [0,1], 'splitting', 'on');
dduf = chebfun(ddu, [0,1], 'splitting', 'on');
ddduf = chebfun(dddu, [0,1], 'splitting', 'on');

Tf = chebfun(T, [0,1], 'splitting', 'on');
dTf = chebfun(dT, [0,1], 'splitting', 'on');
ddTf = chebfun(ddT, [0,1], 'splitting', 'on');

cf = chebfun(c, [0,1], 'splitting', 'on');
dcf = chebfun(dc, [0,1], 'splitting', 'on');
ddcf = chebfun(ddc, [0,1], 'splitting', 'on');

figure(1)
plot(uf)
xlabel('y'), ylabel('u')
title('x velocity')
%print -dpdf out.pdf
print -dpsc out.ps

figure(2)
plot(duf)
xlabel('y'), ylabel('u_y')
title('First derivative of x velocity')
%print -dpdf -append out.pdf
print -dpsc -append out.ps

figure(3)
plot(dduf)
xlabel('y'), ylabel('u_{yy}')
title('Second derivative of x velocity')
%print -dpdf -append out.pdf
print -dpsc -append out.ps

figure(4)
plot(ddduf)
xlabel('y'), ylabel('u_{yyy}')
title('Third derivative of x velocity')
%print -dpdf -append out.pdf
print -dpsc -append out.ps

figure(5)
plot(Tf)
xlabel('y'), ylabel('T')
title('Temperature')
%print -dpdf out.pdf
print -dpsc -append out.ps

figure(6)
plot(dTf)
xlabel('y'), ylabel('T_y')
title('First derivative Temperature')
%print -dpdf out.pdf
print -dpsc -append out.ps

figure(7)
plot(ddTf)
xlabel('y'), ylabel('T_{yy}')
title('Second derivative Temperature')
%print -dpdf out.pdf
print -dpsc -append out.ps

figure(8)
plot(cf)
xlabel('y'), ylabel('s')
title('Concentration')
%print -dpdf out.pdf
print -dpsc -append out.ps

figure(9)
plot(dcf)
xlabel('y'), ylabel('s_y')
title('First derivative concentration')
%print -dpdf out.pdf
print -dpsc -append out.ps

figure(10)
plot(ddcf)
xlabel('y'), ylabel('s_{yy}')
title('Second derivative concentration')
%print -dpdf out.pdf
print -dpsc -append out.ps

save('baseflow.mat','uf','duf','dduf','ddduf',...
                    'Tf','dTf','ddTf',...
                    'cf','dcf','ddcf')
fprintf(1,'Saved chebyshev functions into baseflow.mat\n')
