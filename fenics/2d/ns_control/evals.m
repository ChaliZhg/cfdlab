load state.mat;

%sym = 'o';
A = A - B*S; sym = 'r*';

n = 50;
opts.p = 4*n;

[V,D,flag] = eigs(A,M,n,'SM',opts);
assert(flag==0)
D = diag(D);
plot(real(D),imag(D),sym,'LineWidth',1.5)
grid on
D(1:10)
