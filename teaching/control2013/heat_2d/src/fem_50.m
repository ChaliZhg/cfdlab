function fem_50 ( )
   parameters;
%
%  Read the nodal coordinate data file.
%
  load coordinates.dat;
%
%  Read the triangular element data file.
%
  eval ( 'load elements3.dat;', 'elements3=[];' );
%
%  Read the Neumann boundary condition data file.
%
  eval ( 'load neumann.dat;', 'neumann=[];' );
%
%  Read the Dirichlet boundary condition data file.
%
  eval ( 'load dirichlet.dat;', 'dirichlet=[];' );


  nNodes = size(coordinates,1)
  nTri = size(elements3,1)

% Observations on left boundary
  nObs = 3;
  H = sparse ( nObs, nNodes );
  ymid = 0.5*( coordinates(neumann(:,1),2) + coordinates(neumann(:,2),2) );
  ymin = [0.20, 0.50, 0.80];
  ymax = [0.25, 0.55, 0.85];
  for j=1:nObs
    ie = find(ymid > ymin(j) & ymid < ymax(j));
    ds = abs(coordinates(neumann(ie,2),2) - coordinates(neumann(ie,1),2));
    H(j,neumann(ie,1)) = H(j,neumann(ie,1)) + 0.5*ds';
    H(j,neumann(ie,2)) = H(j,neumann(ie,2)) + 0.5*ds';
    ds = sum(ds);
    H(j,:) = H(j,:)/ds;
  end

  A = sparse ( nNodes, nNodes );
  M = sparse ( nNodes, nNodes );
%
%  Assembly.
%
  for j = 1 : nTri
%   Stiffness matrix
    A(elements3(j,:),elements3(j,:)) = A(elements3(j,:),elements3(j,:)) ...
      + nu*stima3(coordinates(elements3(j,:),:));
%   Mass matrix
    M(elements3(j,:),elements3(j,:)) = M(elements3(j,:), elements3(j,:)) ...
       + det([1,1,1;coordinates(elements3(j,:),:)'])*[2,1,1;1,2,1;1,1,2]/24;
  end
%
%  Determine which nodes are associated with Dirichlet conditions.
%
  BoundNodes = unique ( dirichlet );
  ControlNodes = find ( coordinates(:,1) - a > -eps );
  FreeNodes = setdiff ( 1:nNodes, BoundNodes );

  B = -A(FreeNodes, ControlNodes);
  A = -A(FreeNodes, FreeNodes);
  M = M(FreeNodes, FreeNodes);

  nBoundNodes = size(BoundNodes,1)
  nFreeNodes = size(FreeNodes,2)
  nControlNodes = size(ControlNodes,1)

% Add shift
  A = A + omega * M;
  B = B * sin(pi*coordinates(ControlNodes,2));

  A1 = M \ A;
  B1 = M \ B;
  fprintf(1,'Done inversion\n')
  [X,L,G] = care(A1, B1, M, 1, [], M);
  fprintf(1,'Done care\n')

% Compute eigenvalues
  e = eigs(A, M, 20, 'LR')
  figure(10)
  plot(real(e),imag(e),'o')
  pause

  Z = sparse ( nNodes, 1 );        % Complete solution
  z = sparse ( nFreeNodes, 1 );    % Interior solution
  u = sparse ( nControlNodes, 1 ); % Control vector

  dt = 0.01;
  Nt = 100;
  A1 = M - dt*A;
  Z(:)= sin(0.5*pi*coordinates(:,1)) .* sin(pi*coordinates(:,2));
  z(:) = Z(FreeNodes);
  u(:) = Z(ControlNodes);

% Time loop
  for it = 1:Nt
     z = A1 \ (M*z + dt*B*u);
     Z(FreeNodes) = z;
     Z(ControlNodes) = u;
     show ( elements3, coordinates, full ( Z ) );
     pause(1)
  end
%
%  Graphic representation.
%
  show ( elements3, coordinates, full ( Z ) );

end
