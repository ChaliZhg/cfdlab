function batch_vortex()

global U  k

x0=0;
y0=0;
Re = 800;
a = 1.0; % vortex core radius
U = 1.0;

k = 3.8317/a;
nu = U*a/Re;

xmin=-2*a;
xmax=+2*a;
ymin=-2*a;
ymax=+2*a;

n=100;
x=linspace(xmin,xmax,n);
y=linspace(ymin,ymax,n);
[X,Y]=meshgrid(x,y);
omg = zeros(n,n);
u = zeros(n,n);
v = zeros(n,n);

J1 = @(r) besselj(1, r);

for i=1:n
   for j=1:n
      r = sqrt( (X(i,j)-x0)^2 + (Y(i,j) - y0)^2 );
      theta = atan2(Y(i,j), X(i,j));
      if r <= a
         %omg(i,j) = C * k^2 * X(i,j) * J1(k*r) / r;
         omg(i,j) = U * k^2 * X(i,j) * J1(k*r) / r;
         [u(i,j),v(i,j)] = vel(r, theta);
      else
         omg(i,j) = 0;
         u(i,j) = 0;
         v(i,j) = 0;
      end
   end
end

figure(1)
omin = min(min(omg))
omax = max(max(omg))
contourf(X, Y, omg, linspace(omin,omax,20))
colorbar

figure(2)
contourf(X, Y, u)
colorbar

figure(3)
contourf(X, Y, v)
colorbar

end

%-----------------------------------------------------------------------------
function [u,v] = vel(r, theta)

global U  k

kr = k*r;
f = besselj(0,kr) - besselj(1,kr)/kr;
g = sin(theta)^2 * U * besselj(1,kr)/r;

u = U * sin(theta) * f * k * cos(theta) - g;
v = -cos(theta)^2 * U * f * k - g;

end
