load Mc.mat;
load Ac.mat;

[V,D] = eigs(Ac,Mc,50,'SM');
D = diag(D);
plot(real(D),imag(D),'o','LineWidth',1.5)
grid on
D(1:10)
