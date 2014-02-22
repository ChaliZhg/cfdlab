load linear.mat;

[V,D] = eigs(A,M,50,'SM');
D = diag(D);
plot(real(D),imag(D),'o','LineWidth',1.5)
grid on
D(1:10)
