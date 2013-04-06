clear all
close all

parameters;
[A,B] = get_system_matrices(m,M,k,d,I,g,l,alpha,beta);

% condition initiale assez loin du point d'équilibre instable (position
% verticale)  \theta0 = 0.2 rad soit environ 11.45 degré
x0 = [0.5; 0.2; 0.4; 1];
tspan = [0:0.01:20];
options = odeset('RelTol',1e-8,'AbsTol',1e-8);

[t,x] = ode15s(@fbo,tspan,x0,options,M,m,l,g,k,d,I) ;

% BO , pas de contrôle
figure(1),  title('Evolution of state variables without control'),
         subplot(2,2,1), plot(t,x(:,1)), title('Position of cart'),
         subplot(2,2,2) ; plot(t,x(:,2)) ; title('Speed of cart') ;
         subplot(2,2,3) ; plot(t,x(:,3)) ; title('Angle of pendulum') ;
         subplot(2,2,4), plot(t,x(:,4)), title('Angular speed of pendulum') ;

% BF contrôle en BF         
C = [1 0 0 0; 
     0 0 1 0];
Q = C'*C;
[X,L,K] = care(A, B, Q);
[t,x] = ode15s(@fbf,tspan,x0,options,M,m,l,g,k,d,I,K,alpha,beta) ;

figure(2),  title('Evolution of state variables with control'),
         subplot(2,2,1), plot(t,x(:,1)), title('Position of cart'),
         subplot(2,2,2) ; plot(t,x(:,2)) ; title('Speed of cart') ;
         subplot(2,2,3) ; plot(t,x(:,3)) ; title('Angle of pendulum') ;
         subplot(2,2,4), plot(t,x(:,4)), title('Angular speed of pendulum') ;
