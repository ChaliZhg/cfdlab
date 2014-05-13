clear all
load linear.mat
load freeinds.txt
load pinds.txt
who

nc = size(B,2);   % number of control variables
nu = 4;           % how many eigenvalues to compute
shift = 0;

% number of Lanczos vectors
opts.p = 50;

[Vt,D1] = eigs(A,M,nu,'SM',opts);
disp('Eigenvalues of A')
D1=diag(D1)

[Zt,D2] = eigs(A',M',nu,'SM',opts);
disp('Eigenvalues of A^T')
D2=diag(D2)

% find unstable eig
iu = find(real(D1) > 0);
nu = length(iu);
fprintf(1, 'Number of unstable eigenvalues of A = %d\n', nu)

D1 = D1(iu)
Vt = Vt(:,iu);

% find unstable eig
iu = find(real(D2) > 0);
nu = length(iu);
fprintf(1, 'Number of unstable eigenvalues of A^T= %d\n', nu)

D2 = D2(iu)
Zt = Zt(:,iu);

disp('Is the order of eigenvalues same ? Enter to continue')
pause

% NOTE: check that eigenvalues are in same order

% make Vt and Zt orthonormal
% p must be diagonal
disp('Following must be a diagonal matrix. Is it ?')
p = Vt.' * M * Zt
p = diag(p);

% normalize
for j=1:nu
   Zt(:,j) = Zt(:,j) / p(j);
end

% freeinds, pinds are indices inside fenics
% We have to shift by one since python indexing starts at 0 but matlab 
% starts at 1
freeinds = freeinds + 1;
pinds    = pinds + 1;
% get matlab indices of velocity+temperature
[tmp,vTinds] = setdiff(freeinds, pinds, 'stable');
% get matlab indices of pressure
pinds = setdiff(1:length(freeinds), vTinds, 'stable');

% eigenvector component for velocity+temperature
Vty = Vt(vTinds,:);  % eigenvectors of (A,M)
Zty = Zt(vTinds,:);  % eigenvectors of (A',M')
Ztp = Zt(pinds,:);

E11 = M(vTinds,vTinds);
A11 = A(vTinds,vTinds);
A12 = A(vTinds,pinds);
B1  = B(vTinds,:);
B2  =-B(pinds,:);

% check orthonormality
disp('Is this identity matrix ?')
p = Vty.' * E11 * Zty

U = (1/sqrt(2)) * [1,   1; ...
                   1i, -1i];

Vy = Vty * U';
Zy = Zty * U.';
Zp = Ztp * U.';


disp('Vy and Zy must be real')
max(abs(imag(Vy)))
max(abs(imag(Vy)))

% Vy and Zy must be real, making sure imaginary part is close to zero
Vy = real(Vy);
Zy = real(Zy);
Zp = real(Zp);

disp('Is this identity matrix ?')
p = Vy.' * E11 * Zy


% Compute B12
np = length(pinds);
ny = length(vTinds);
N  = [E11, A12; A12' sparse(np,np)];
RHS= [sparse(ny,nc); B2];
Z1 = N\RHS;
B12= B1 + A11*Z1(1:ny,:);

% Project to unstable subspace
Au = Zy' * A11 * Vy;
Bu = Zy' * B12;
Ru = eye(nc);

% minimal norm control
%Qu = zeros(size(Au));

% LQR problem
Qu = Vy' * E11 * Vy;

[Pu,L,G]=care(Au,Bu,Qu,Ru);
disp('Eigenvalues of projected system with feedback')
L

B = sparse([B1; -B2]);
E11 = sparse(E11);
Z = sparse([Zy; Zp]);
Zy = sparse(Zy);
Pu = sparse(Pu);
Kt = (B' * Z) * Pu * (Zy' * E11);
S  = [Kt, sparse(nc,np)];
A  = sparse([A11, A12; A12', sparse(np,np)]);
M  = sparse([E11, sparse(ny,np); sparse(np,ny+np)]);
[V,D] = eigs(A-B*S,M,nu,'SM',opts);
disp('Eigenvalues of full system with feedback')
diag(D)

Kt = full(Kt);
save('gain.mat','Kt')

save('state.mat','M','A','B','S','Z')
