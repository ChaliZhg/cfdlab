% p = degree of bspline
% N = number of cells = no. of knot intervals
% cflmode = 'zhang' or 'praveen'

function [Ul, ubar, mintheta] = project(U)

globals;

[N,p] = size(U);
p = p - 1;

mintheta = 1.0e20;

% limiter
for j=1:N
   ubar(j) = sum(U(j,:))/(p+1);
   ff = bezier(U(j,:), xx);
   mj = min(ff);
   Mj = max(ff);
   a1 = abs(umax - ubar(j))/(abs(Mj - ubar(j)) + 1.0e-14);
   a2 = abs(umin - ubar(j))/(abs(mj - ubar(j)) + 1.0e-14);
   theta = min( [a1, a2, 1.0] );
   Ul(j,:) = theta*(U(j,:) - ubar(j)) + ubar(j);
   mintheta = min([mintheta, theta]);
end
