load linear.mat

[V1,D1] = eigs(A,M,2,'SM');
diag(D1)

[V2,D2] = eigs(A',M',2,'SM');
diag(D2)

B = [Bv, Bt, Bh];


% check controllability
B' * V2(:,1)
B' * V2(:,2)

% number of unstable eigenvalues
nu = 2;
Du = D1(1:nu,1:nu);
