clear all
close all

parameters;
[A,B] = get_system_matrices(m,M,k,d,I,g,l,alpha,beta);

% Eigenvalues of un-controlled system
eo = eig(A);
plot(real(eo), imag(eo), 'o')
figure(1)
grid on

% minimal norm control
% we need to shift eigenvalues of A since there is one zero eigenvalue
om=0.1;
[X,L,K] = care(A+om*eye(size(A)), B, zeros(size(A)), 1);
L=eig(A-B*K);
figure(2)
plot(real(eo), imag(eo), 'o', real(L), imag(L), '*')
grid on

% Stabilize all using LQR
C = eye(4); 
Q = C'*C;
[X,L,K] = care(A, B, Q, 1);
figure(3)
plot(real(eo), imag(eo), 'o', real(L), imag(L), '*')
grid on

% Stabilize only position and angle using LQR
C = [1 0 0 0; 
     0 0 1 0];
Q = C'*C;
[X,L,K] = care(A, B, Q, 1);
figure(4)
plot(real(eo), imag(eo), 'o', real(L), imag(L), '*')
grid on
