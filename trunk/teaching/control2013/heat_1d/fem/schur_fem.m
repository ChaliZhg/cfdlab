function K = schur_fem(M,A,B,C,R,D,N)
% méthode de Schur pour résoudre l'ARE avec schéma EF

npts = size(B,1);
n = npts+1;
h = 1/n;

% avec lqr
[K,X] = lqr(M\A, M\B, C'*C, R+D^2, N);


K = real(K)/M;
% noyau sur [0,1]
k_EF = [0 K 0];


x = 0:h:1;
figure(1),
plot(x,k_EF,'o-r'),hold on
xlabel('x','fontsize',24),ylabel('k(x)','fontsize',24),
legend('noyau k en fct de x')

end

