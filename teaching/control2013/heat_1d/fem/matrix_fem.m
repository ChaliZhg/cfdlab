function [M,A,B,Q]=matrix_fem(n,nu,alpha)
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

% mass matrice
M = h*spdiags([e/6 2*e/3 e/6], -1:1, n-1, n-1);

A = nu*spdiags([e -2*e e], -1:1, n-1, n-1)/h + alpha*M;

% matrice de sortie
Q = M;

B = sparse(n-1,1);
B(n-1,1) = nu/h;
