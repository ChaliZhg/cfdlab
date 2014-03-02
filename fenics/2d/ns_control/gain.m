clear all
load linear.mat

nu = 2;

% number of Lanczos vectors
opts.p = 50;

[V1,D1] = eigs(A,M,nu,'SM',opts);
disp('Eigenvalues of A')
diag(D1)

[V2,D2] = eigs(A',M',nu,'SM',opts);
disp('Eigenvalues of A^T')
diag(D2)

% NOTE: check that eigenvalues are in same order

% make V1 and V2 orthonormal
% p must be diagonal
disp('Is this diagonal')
p = conj(V2')*M*V1
p = diag(p);

% normalize
for j=1:nu
   V2(:,j) = V2(:,j) / p(j);
end
% check orthonormality
disp('Is this diagonal')
p = conj(V2')*M*V1

B = [Bv, Bt, Bh];


% check controllability by hautus
for j=1:nu
   B' * V2(:,j)
end

% number of unstable eigenvalues
Du = D1(1:nu,1:nu);
