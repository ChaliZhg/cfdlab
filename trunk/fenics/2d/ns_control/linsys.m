clear all
load state.mat;

energy = @(x) full(0.5*x'*M*x);

% Intial condition using unstable eigenvalue
x1 = 0.01 * Z(:,1);

fprintf(1,'Initial energy = %e\n', energy(x1))

N = length(x1);

% time step
dt = 0.1;
iter = 0;
t = 0;
Tf = 200;

% first step: BDF1
% M(x(n) - x(n-1))/dt = A x(n) - B K x(n-1)
% x = x(n), x1 = x(n-1)
M1 = M/dt - A;
A1 = M/dt - B*S;
x = M1 \ (A1 * x1);
x2 = x1;
x1 = x;
iter = iter + 1;
t(iter) = dt;
e(iter) = energy(x);

fprintf(1,'%d %e %e\n',iter,t(iter),e(iter))

% second step onwards: BDF2
% M(1.5*x(n) - 2*x(n-1) + 0.5*x(n-2))/dt = A x(n) - B K x(n-1)
% x = x(n), x1 = x(n-1), x2 = x(n-2)
M1 = 1.5*M/dt - A;
A1 = 2*M/dt - B*S;
A2 =-0.5*M/dt;
[L1,U1,P1,Q1] = lu(M1);

while t<Tf
   rhs = A1*x1 + A2*x2;
   x = Q1 * (U1 \ (L1 \ (P1 * rhs ) ) );
   x2 = x1;
   x1 = x;
   iter = iter + 1;
   t(iter) = t(iter-1) + dt;
   e(iter) = energy(x);
   fprintf(1,'%d %e %e\n',iter,t(iter),e(iter))
   if(mod(iter,100)==0) 
      semilogy(t,e)
      pause(0.1)
   end
end
