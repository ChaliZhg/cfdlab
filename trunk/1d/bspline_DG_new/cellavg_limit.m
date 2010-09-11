% p = degree of bspline
% N = number of cells = no. of knot intervals
% cflmode = 'zhang' or 'praveen'

function [Ul, ubar, mintheta] = cellavg_limit(U)

globals;

[N,p] = size(U);
p = p - 1;

mintheta = 1.0e20;

for j=1:N
   ubar(j) = sum(U(j,:))/(p+1);
end

% Find min-max for each cell
%minmax(U,ubar);

% limiter
for j=1:N

   % Theta limiter of Zhang-Shu
   ff = bezier(U(j,:), xx);
   mj = min(ff);
   Mj = max(ff);
   a1 = abs(umax(j) - ubar(j))/(abs(Mj - ubar(j)) + 1.0e-14);
   a2 = abs(umin(j) - ubar(j))/(abs(mj - ubar(j)) + 1.0e-14);
   theta = min( [a1, a2, 1.0] );
   mintheta = min([mintheta, theta]);
   Ul(j,:) = theta*(U(j,:) - ubar(j)) + ubar(j);

end
