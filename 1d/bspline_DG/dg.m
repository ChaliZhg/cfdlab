%clear all
%close all

function [ndof, err] = dg(p,N)

globals;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set test case
%testcase = burger_step;
%testcase = burger_sine;
testcase = lincon_sine


% degree of bspline
%p=1

% number of cells = no. of knot intervals
%N = 50
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Control variables
U    = zeros(N,p+1);
Uold = zeros(N,p+1);
res  = zeros(N,p+1);

% burgers with initial shock -- step function
if testcase==burger_step
   fluxfun = burger;
   % domain size
   xmin =-1.0
   xmax = 1.0
   tfinal = 1.0;

   periodic = no
   for j=1:N/2
      U(j,:) = 1.0;
   end
elseif testcase==burger_sine
   fluxfun = burger;
   % domain size
   xmin =-1.0
   xmax = 1.0
   tfinal = 2.0;

   periodic = yes;
elseif testcase==lincon_sine
   fluxfun = lincon;
   % domain size
   xmin = 0.0
   xmax = 1.0
   tfinal = 2.0;
   periodic = yes;
end

% Cell size
h = (xmax - xmin)/N;

% weights and quadrature points for limiting - zhang-shu
if p==1
   xx = [-0.5 +0.5];
   cfl= 1.0/4.0;
elseif p==2 || p==3
   xx = [-0.5 0.0 +0.5];
   cfl= 1.0/6.0;
elseif p==4 || p==5
   xx = [-0.5 -1/sqrt(20) +1/sqrt(20) +0.5];
   cfl= 1.0/12.0;
end
xx = xx + 0.5; % shift to [0,1]

% cell centers
for j=1:N
   xc(j) = xmin + (j-1)*h + 0.5*h;
end

% mass matrix
for m=0:p
   for n=0:p
      M(m+1,n+1) = nchoosek(p,m) * nchoosek(p,n) * ...
                   gamma(2*p + 1 - (m+n)) * gamma(m+n+1) / gamma(2+2*p);
   end
end
invM = inv(M);

% initial condition
if testcase ~= burger_step
   for j=1:N
      x1 = xmin + (j-1)*h;
      U(j,1) = initcond(x1);
      U(j,p+1) = initcond(x1+h);
   end

   if p>1
      AA = M(2:p,2:p); invAA = inv(AA);
      b = zeros(p-1,1);
      for j=1:N
         x1 = xmin + (j-1)*h;
         for n=1:p-1
            fun = @(x)( initcond(x1+h*x) .* bernstein(p,n,x) );
            b(n) = quadgk(fun,0,1) - U(j,1)*M(1,n+1) - U(j,p+1)*M(p+1,n+1);
         end
         U(j,2:p) = invAA * b;
      end
   end
end

umin = 1.0e20;
umax =-1.0e20;

% plot initial condition
figure(1)
hold on
for j=1:N
   x1 = xmin + (j-1)*h;
   s  = linspace(0,1,10);
   up = bezier(U(j,:), s);
   xp = x1 + h*s;
   plot(xp,up)
   umin = min(umin, min(up));
   umax = max(umax, max(up));
end
title('Initial condition');

umin
umax

% Parameters for runge-kutta
ark = [0.0 3.0/4.0 1.0/3.0];
brk = 1 - ark;

amax  = 1; % max speed
dt    = cfl*h/amax;
niter = tfinal/dt;

% Cell average values
for j=1:N
   ubar(j) = sum(U(j,:))/(p+1);
end

time = 0.0;

% time integration
for iter=1:niter
   Uold = U;

   % Runge-kutta stages
   for rks=1:3

      Ul       = U;
      res(:,:) = 0.0;

      % limiter
      for j=1:N
         ff = bezier(U(j,:), xx);
         mj = min(ff);
         Mj = max(ff);
         a1 = abs(umax - ubar(j))/abs(Mj - ubar(j));
         a2 = abs(umin - ubar(j))/abs(mj - ubar(j));
         theta = min( [a1, a2, 1.0] );
         Ul(j,:) = theta*(U(j,:) - ubar(j)) + ubar(j);
      end


      % Inter-element fluxes
      if periodic==no
         flx = numflux(Ul(1,1), Ul(1,1));
      else
         flx = numflux(Ul(N,p+1), Ul(1,1));
      end
      res(1,1) = res(1,1) - flx;

      for j=1:N-1
         flx = numflux(Ul(j,p+1), Ul(j+1,1));
         res(j,  p+1) = res(j,  p+1) + flx;
         res(j+1,1  ) = res(j+1,1  ) - flx;
      end

      if periodic==no
         flx = numflux(Ul(N,p+1), Ul(N,p+1));
      else
         flx = numflux(Ul(N,p+1), Ul(1,1));
      end
      res(N,p+1) = res(N,p+1) + flx;

      % Flux integral
      for j=1:N
         for n=0:p
            fun = @(x)( flux( bezier(Ul(j,:), x) ).*dbernstein(p,n,x) );
            res(j,n+1) = res(j,n+1) - quadgk(fun,0,1);
         end
      end

      % Update control variables
      for j=1:N
         res(j,:) = invM * (res(j,:))';
         U(j,:) = ark(rks)*Uold(j,:) + ...
                  brk(rks)*(Ul(j,:) - (dt/h)*res(j,:));
      end

   end % end of RK


   % Cell average values
   for j=1:N
      ubar(j) = sum(U(j,:))/(p+1);
   end

   % limiter
   mintheta=1.0e20;
   for j=1:N
      ff = bezier(U(j,:), xx);
      mj = min(ff);
      Mj = max(ff);
      a1 = abs(umax - ubar(j))/abs(Mj - ubar(j));
      a2 = abs(umin - ubar(j))/abs(mj - ubar(j));
      theta = min( [a1, a2, 1.0] );
      mintheta = min(mintheta,theta);
      Ul(j,:) = theta*(U(j,:) - ubar(j)) + ubar(j);
   end

   time = time + dt;
   fprintf(1,'%d %f %f %f\n', iter, dt, time, mintheta);

   figure(2)
   for j=1:N
      x1 = xmin + (j-1)*h;
      s  = linspace(0,1,10);
      up = bezier(Ul(j,:), s);
      xp = x1 + h*s;
      plot(xp,up)
      hold on
   end
   grid
   plot(xc,ubar,'o')
   hold off

end


% Compute error for linear convection equation
err = 0.0;
for j=1:N
   x1  = xmin + (j-1)*h;
   fun = @(x)( (bezier(U(j,:), x) - initcond(x1+h*x-time)).^2 );
   err1= quadgk(fun,0,1,'AbsTol',1e-10,'RelTol',0);
   err = err + err1;
end
err = sqrt(h*err);

fprintf(1,'L2 error = %f\n', err);

ndof = N*(p+1);
