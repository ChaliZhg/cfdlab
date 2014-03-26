clear all
load linear.mat
load freeinds.txt
load pinds.txt
who

nu = 2;
shift = 0;

% number of Lanczos vectors
opts.p = 50;

[V1,D1] = eigs(A,M,nu,'SM',opts);
disp('Eigenvalues of A')
diag(D1)

[V2,D2] = eigs(A',M',nu,'SM',opts);
disp('Eigenvalues of A^T')
diag(D2)

% NOTE: check that eigenvalues are in same order

% freeinds, pinds are indices inside fenics
% We have to shift by one since python indexing starts at 0 but matlab 
% starts at 1
freeinds = freeinds + 1;
pinds    = pinds + 1;
% get matlab indices of velocity+temperature
[tmp,vTinds] = setdiff(freeinds, pinds);
% get matlab indices of pressure
pinds = setdiff(1:length(freeinds), vTinds);

% eigenvector component for velocity+temperature
V1y = V1(vTinds,:);  % eigenvectors of (A,M)
V2y = V2(vTinds,:);  % eigenvectors of (A',M')

E11 = M(vTinds,vTinds);
A11 = A(vTinds,vTinds);
A12 = A(vTinds,pinds);
B1  = B(vTinds,3);
B2  = B(pinds,3);

% make V1 and V2 orthonormal
% p must be diagonal
disp('Is this diagonal matrix ?')
p = conj(V2')*M*V1
p = diag(p);

% normalize
for j=1:nu
   V2(:,j) = V2(:,j) / p(j);
end
% check orthonormality
disp('Is this identity matrix ?')
p = conj(V2')*M*V1

% check controllability by hautus
for j=1:nu
   B.' * V2(:,j)
end

% number of unstable eigenvalues
Du = D1(1:nu,1:nu) + shift*eye(nu);
% Stabilization
Bu=V2.'*B;
Ru=eye(size(Bu,2));
Qu=zeros(nu);
[Pu,L,G]=care(Du,Bu,Qu,Ru);
disp('eigenvalues with feedback of projected system')
L
K=(((B.')*V2)*Pu)*(V1.'*M);
disp('norm of imag(K)')
norm(imag(K))
K=real(K);
K=sparse(K);
A=sparse(A);
M=sparse(M);
B=sparse(B);
[Vc,Dc]=eigs(A-B*K,M,5,'SM',opts);
disp('eigenvalues with feedback of full system')
diag(Dc)
