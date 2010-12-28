% Program to find and plot convergence rate of CC and RE
% Ignore first iskip data points
function rate = conv_rate (iskip)

if nargin==0
   istart = 1;
else
   istart = iskip + 1;
end

data = load('error.dat');

[ns,nf] = size(data);
nf = nf - 2;

x    = data(istart:ns, 2);
xall = data(:,2);
x1   = xall(1);
x2   = xall(ns);
for j=1:nf
   f    = data(istart:ns, j+2);
   fall = data(:, j+2);
   A = [ones(size(x)), log10(x)];
   sol = inv(A'*A)*A'*log10(f);
   rate(j) = sol(2);
   f1(j) = 10^(sol(1) + sol(2)*log10(x1));
   f2(j) = 10^(sol(1) + sol(2)*log10(x2));
   if j==1
      loglog(xall,fall,'o','MarkerSize',8)
   else
      loglog(xall,fall,'*','MarkerSize',8)
   end
   hold on
end
legend('CC','RE')
for j=1:nf
   loglog([x1 x2],[f1(j) f2(j)],'-','LineWidth',1.5)
end
hold off
