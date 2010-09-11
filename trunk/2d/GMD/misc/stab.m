n  = 100;
lx = linspace(-pi,+pi,n);
ly = linspace(-pi,+pi,n);
[X,Y]=meshgrid(lx,ly);

for nux=0.1:0.1:1.1
   Z = 1 - 0.5*nux*(1 + cos(Y)).*(1 - exp(-i*X));
   Z = abs(Z);
   fprintf(1,'%e %e\n', nux, max(max(Z)));
   contour(X,Y,Z);
   pause
end

