function fem_50 ( )
   globals;
   parameters;

   read_data();

   [M,A,B,Q,H] = get_matrices();

   % Compute feedback gain matrix

   % Compute estimation gain matrix

% Compute eigenvalues and eigenfunctions
% [V,D] = eigs(A, M, 5, 'LR');
% e = diag(D);
% figure(10)
% plot(real(e),imag(e),'o')
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
  Nt = 100;
  A1 = M - dt*A;
  za(:)= -cos(0.5*pi*coordinates(:,1)) .* sin(pi*coordinates(:,2));
  z(:) = za(FreeNodes);
  u(:) = za(ControlNodes);
  v    = 0;
  show ( elements3, coordinates, full(za) );
  pause(2)

% Time loop
  t = 0;
  energy = zeros(Nt,1);
  time = zeros(Nt,1);
  for it = 1:Nt
     z = A1 \ (M*z + dt*B*v);
     t = t + dt;
     za(FreeNodes) = z;
     za(ControlNodes) = v*u;
     time(it) = t;
     energy(it) = z'*M*z;
     y = H*z;
     fprintf(1,'Time = %e, Energy = %e, Obs = %e %e %e\n', ...
             t, energy(it), y(1), y(2), y(3))
     if mod(it,10) == 0
      show ( elements3, coordinates, full(za) );
      pause(1)
     end
  end
  show ( elements3, coordinates, full(za) );

  % Plot energy
  figure(10)
  semilogy(time, energy)
  xlabel('Time')
  ylabel('Energy')

end
