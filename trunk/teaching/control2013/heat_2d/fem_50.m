function fem_50 ( )
   globals;
   parameters;

   read_data();

   [M,A,B,Q,H] = get_matrices();
   
%==========================================================================
% EVALUATION OF THE FEEDBACK MATRIX
%==========================================================================
   
%--------------------------------------------------------------------------
% Feedback matrix based on the whole system
%--------------------------------------------------------------------------
%      R = 1;
%      [X,L,K] = care(full(A),full(B), full(Q),R,[],full(M));
%      A=A-B*sparse(K); 

%--------------------------------------------------------------------------
% Feedback matrix based on only the unstable components
%--------------------------------------------------------------------------
   % Unstable eigenvalues and eigenvectors
   nu = 1;
   [V,D] = eigs(A,M,nu,'lr');
   % Bu is the matrix acting on the control corresponding to the unstable 
   % portion of the diagonalized system, 
   Bu = V'*B;

   % Solving the ARE for the unstable part of the diagonalized system
   Qu = zeros(nu);
   Ru = 1;
   [Pu,L,G] = care(D,Bu,Qu,Ru);

   disp('The initial unstable eigenvalues were')
   eig(D)
   disp('The unstable eigenvalues are modified to')
   eig(D - Bu*Bu'*Pu)

   % Matrix P for the initial system.
   P = sparse(V*Pu*V');

   % Feedback matrix for the original system
   K = B'*P*M;
%--------------------------------------------------------------------------   

% Compute estimation gain matrix

% Compute eigenvalues and eigenfunctions
 %[V,D] = eigs(A-B*K, M, 5, 'LR');
 %e = diag(D);
 %figure(11)
 %plot(real(e),imag(e),'o')
% fprintf(1,'Most unstable eigenvalue (numerical) = %f\n', e(1))
% fprintf(1,'Most unstable eigenvalue (exact)     = %f\n', -pi^2/40 + omega)
% ef = zeros(nNodes,1);
% ef(FreeNodes) = V(:,1);
% show ( elements3, coordinates, full ( ef ) );
% pause(2)

  za = sparse ( nNodes, 1 );        % Complete solution
  z  = sparse ( nFreeNodes, 1 );    % Interior solution
  u  = sparse ( nControlNodes, 1 ); % Control vector

  dt = 0.01;
  Nt = 10000;
  A1 = M - dt*(A-B*K);
  za(:)= -cos(0.5*pi*coordinates(:,1)) .* sin(pi*coordinates(:,2));
  z(:) = za(FreeNodes);
  u(:) = sin(pi*coordinates(ControlNodes,2)); % sin(pi*y)
  show ( elements3, coordinates, full(za) );
  pause(2)

% Time loop
  t = 0;
  energy = zeros(Nt,1);
  time = zeros(Nt,1);
  for it = 1:Nt
     z = A1 \ (M*z);
     t = t + dt;
     za(FreeNodes) = z;
     v = -K*z;
     za(ControlNodes) = v*u;
     time(it) = t;
     energy(it) = z'*M*z;
     y = H*z;
     fprintf(1,'Time = %e, Energy = %e, Obs = %e %e %e\n', ...
             t, energy(it), y(1), y(2), y(3))
     if mod(it,500) == 0
      show ( elements3, coordinates, full(za) );
      pause(1)
     end
     v = -K*z;
  end
  show ( elements3, coordinates, full(za) );

  % Plot energy
  figure(10)
  semilogy(time, energy)
  xlabel('Time')
  ylabel('Energy')

end
