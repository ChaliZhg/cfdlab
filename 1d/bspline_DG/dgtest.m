clear all
close all

globals;
init();

N = [20 40 80 160];
p    = 1;

ndof = [];
err  = [];
for j=1:length(N)
   [a1, a2] = dg(p,N(j));
   ndof = [ndof a1];
   err  = [err a2];
   figure(10)
   loglog(ndof,err,'o-');
end

for j=2:length(ndof)
   pp =-(log10(err(j)) - log10(err(j-1)))/(log10(ndof(j)) - log10(ndof(j-1)))
end
