clear all
close all

parameters;
[A,B] = get_system_matrices(m,M,k,d,I,g,l,alpha,beta);

Q = diag([1/0.68^2, 1, 1/0.175^2, 1/0.5^2]);
Ru = 1/3^2;

[K,X] = lqr(A, B, Q, Ru);
disp('Eigenvalues of A-B*K')
eig(A-B*K)

% Observation operator
H = [1, 0, 0, 0;
     0, 0, 1, 0];

% Number of time steps
Nt = 1500;

% Initial condition
x0 = [0.5; 0.2; 0.4; 1];

% Noise in state
w = (5e-1)^2 * randn(Nt,4);
Rw = diag(std(w).^2);

% Noise in observation
v = (5e-2)^2 * randn(Nt,2);
Rv = diag(std(v).^2);

[L,S] = lqr(A', H', Rw, Rv);
L = real(L');

disp('Eigenvalues of A-L*H')
eig(A-L*H)

% Coupled system
Ae = [A,   -B*K; ...
      L*H,  A-L*H-B*K];

disp('Eigenvalues of coupled system')
eig(Ae)

% Noise matrix
Be = [eye(4,4),   zeros(4,2); ...
      zeros(4,4), L];

% Observation matrix
He = [H,          zeros(2,4); ...
      zeros(2,4), H];

% N : number of state variables, no : number of observation
N = 4; no=2;
dt = 0.01;
t = 0:dt:dt*(Nt-1);

znb = zeros(2*N,Nt);
zbc = zeros(2*N,Nt);

u = zeros(1,Nt);

ynb = zeros(4,Nt);
ybc = zeros(4,Nt);

% Initial condition; for estimator we give wrong IC
znb(:,1) = [x0;x0*0];
zbc(:,1) = [x0;x0*0];

% Euler implicite
[L1,U1] = lu(eye(8)-dt*Ae);


for nt=1:Nt-1    
    % calculate state with control, without noise
    % error in initial condition
    %--------------------------------------------
    znb(:,nt+1) = U1\(L1\(znb(:,nt)));
    
    % calculate state with noise
    %---------------------------------------
    zbc(:,nt+1) = U1\(L1\(znb(:,nt) + dt*Be*[w(nt+1,:)'; v(nt+1,:)'] ));
    
    % calculate control
    %--------------------------------------------------------
    u(1,nt+1) = -K*zbc(5:8,nt+1);
    
    % observation
    ynb(:,nt+1) = He*znb(:,nt+1);
    ybc(:,nt+1) = He*zbc(:,nt+1) + [v(nt+1,:)';0;0];
end;

% Validation de l'estimateur en partant d'une mauvaise condition initiale
figure(4),
plot(t,ynb(1,:),t,ynb(3,:),'r')
legend('sortie  correspondant la position du chariot','sortie estimée correspondant à la position du chariot')

% Validation de l'estimateur en partant d'une mauvaise condition initiale
% plus comparaison de l'état bruité et de l'état reconstruit
figure(5)
plot(t,ybc(1,:),t,ybc(3,:),'r')
legend('sortie bruitée correspondant à la position du chariot','sortie estimée correspondant à la position du chariot')

figure(51),
plot(t,zbc(2,:),t,zbc(6,:),'r')
legend('sortie bruitée correspondant à la vitesse du chariot','sortie estimée correspondant à la vitesse du chariot')

figure(6),
plot(t,u),legend('Evolution of control')
