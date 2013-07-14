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

  Z = sparse ( nNodes, 1 );        % Complete solution
  z = sparse ( nFreeNodes, 1 );    % Interior solution
  u = sparse ( nControlNodes, 1 ); % Control vector

  dt = 0.01;
  Nt = 100;
  A1 = M - dt*A;
  Z(:)= -cos(pi*coordinates(:,1)) .* sin(pi*coordinates(:,2));
  z(:) = Z(FreeNodes);
  u(:) = Z(ControlNodes);
  v    = 1;
  show ( elements3, coordinates, full ( Z ) );
  pause(2)

% Time loop
  t = 0;
  for it = 1:Nt
     z = A1 \ (M*z + dt*B*v);
     t = t + dt;
     Z(FreeNodes) = z;
     Z(ControlNodes) = v*u;
     if mod(it,10) == 0
      show ( elements3, coordinates, full ( Z ) );
      pause(1)
     end
     y = H*z;
     fprintf(1,'Time = %e %e %e %e\n', t, y(1), y(2), y(3))
  end
%
%  Graphic representation.
%
  show ( elements3, coordinates, full ( Z ) );

end
