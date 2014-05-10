load linear.mat;

n = 50;
opts.p = 2*n;

[V,D] = eigs(A,M,n,'SM',opts);
D = diag(D);
plot(real(D),imag(D),'o','LineWidth',1.5)
grid on
D(1:10)
