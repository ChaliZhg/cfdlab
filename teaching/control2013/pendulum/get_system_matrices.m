function [A,B] = get_system_matrices(m,M,k,d,I,g,l,alpha,beta)
% Stabilisation des eqs du pendule invers�. Eqs d'etat sous la forme
%
%     E x_t = A x + B u + G w
%        y  = H x + Dw w + Du u + v
%
% u : commande, w bruit d'�tat, v bruit de mesure.
%
% INPUT : m : masse du  pendule en kg
%         M : masse du chariot en kg
%         k : coef de frottement du chariot sur le sol
%         d : coef de r�sistance de l'air � la rotation de l'axe
%         I: moment d'inertie kg/m^2
%         l : demi-longueur de l'axe du pendule en m�tres
%         g : intensit� de pesanteur en m/s^2
%     (alpha, beta) permet d'axprimer la force appliqu�e en fonction du
%     voltage V et de dx/dt par la relation : 
%      F(t) = alpha*V(t) - beta*x'(t)
%
% OUTPUT : Matrices d'�tat du pendule invers� associ�es aux vecteurs d'�tat
% 
%         x1 : position of cart
%         x2 : speed of cart
%         x3 : angle du pendule par rapport � la verticale
%         x4 : vitesse angulaire du pendule
%
% Les sorties mesur�es pour l'estimation sont y1 = x1 et y2 = x3.
%
% param�tres � utiliser
%
% M=2.4; m=0.23; g=9.81; l=0.36; I=1/3*m*l^2=0.0099; 
% k = 0.05; d = 0.005;
% alpha = 1.7189,  beta = 7.682
% prendre alpha=1 et beta=0 pour imposer F = force comme contr�le

%% Matrices et syst�me d'�tat

% variables auxilliaires du mod�le 
v1 = (M+m)/(I*(M+m)+l^2*m*M);
v2 = (I+l^2*m)/(I*(M+m)+l^2*m*M);

A = [0   1                      0                       0; ...
     0 -k*v2           -(l*m)^2*g*v2/(I+l^2*m)  l*m*d*v2/(I+l^2*m);...
     0   0                      0                       1; ...
     0 l*m*k*v1/(M+m)        l*m*g*v1                 -d*v1];

% matrice de contr�le
B = [0; v2 ; 0 ;  -l*m*v1/(M+m)];

% CHECK: B contains alpha which should not be added to A
A = A + B*[0 -beta 0 0]

B = B*alpha;
