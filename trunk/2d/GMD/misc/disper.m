N = 100;

%cfl number
nux = 0.8;

kx = linspace(0,2*pi,N);
ky = linspace(0,2*pi,N);

[KX,KY]=meshgrid(kx,ky);
A1 = (1 - 2 * nux * cos(0.5*KY).^2 .* sin(0.5*KX).^2 );
B1 = -nux * cos(0.5*KY).^2 .* sin(KX);
A  = A1 ./ abs(A1 + i*B1);
B  = B1 ./ abs(A1 + i*B1);
disperr = -1./(nux*KX) .* atan(B./A) ;
contourf(KX,KY,disperr)
xlabel('k_1 \Delta x')
ylabel('k_2 \Delta y')
axis equal
axis tight
colorbar

figure(2)
contourf(KX,KY,abs(A1+i*B1))
colorbar

A1 = (1 - 2 * nux * sin(0.5*kx).^2);
B1 = -nux * sin(kx) ;
A  = A1 ./ abs(A1 + i*B1);
B  = B1 ./ abs(A1 + i*B1);
disperr = -1./(nux*kx) .* atan(B./A) ;
figure(3)
plot(kx, disperr)
