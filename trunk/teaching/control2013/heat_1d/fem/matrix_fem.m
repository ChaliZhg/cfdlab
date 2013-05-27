function [M,A,B,C,D,N]=matrix_fem(n,nu,alpha)
% INPUT : n : nombre d'intevalles.
%         nu : viscosité
% Schéma EF P1.
% OUTPUT : Matrices de l'eq de la chaleur sur [0,1] avec CL de Dirichlet,
%          contrôle frontière en x=0.
%
%          M: matrice de masse (identité)
%          A: matrice de rigidité

% pas de la subdivision
h = 1/n;

e = ones(n-1,1);

% matrice de masse
M = h*spdiags([e/6 2*e/3 e/6], -1:1, n-1, n-1);

A = nu*spdiags([e -2*e e], -1:1, n-1, n-1)/h + alpha*M;

% matrice de sortie
C = h*ones(1,n-1);

% matrice de couplage entre état et contrôle.
% D = C(1,1)=h/2;
D = h/2;

N = D*C';

B = sparse(n-1,1);
B(n-1,1) = nu/h;

%R = 1;
