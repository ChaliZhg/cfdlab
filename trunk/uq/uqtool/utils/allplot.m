% Program to find and plot convergence rate of CC and RE
% Ignore first iskip data points
function allplot (istart)

if nargin==0
   istart=1;
end

Jexact = load('Jexact.dat')

Jdata= load('RESULT/J.dat');
ccre = load('RESULT/error.dat');

Ns= Jdata(:,1);
J = Jdata(:,3);
Jc= Jdata(:,4);
h = ccre(:,2);
cc= ccre(:,3);
re= ccre(:,4);

J_error = abs(J - Jexact);
Jc_error= abs(Jc - Jexact);

figure(1)
plot(Ns,J,'o-',Ns,Jc,'o--','LineWidth',1.5)
xlabel('Number of samples')
ylabel('Mean functional')
legend('J','J+CC')

figure(2)
loglog(h,J_error,'o--',...
       h,Jc_error,'*--',...
       h,cc,...
       h,re,'LineWidth',1.5);
xlabel('\Delta\xi')
legend('Error in J','Error in J+CC','CC','RE')

% Convergence rates
ns = length(h);
x = h(istart:ns);

% for J
f = J_error(istart:ns);
A = [ones(size(x)), log10(x)];
sol = inv(A'*A)*A'*log10(f);
rate = sol(2)

% for J+CC
f = Jc_error(istart:ns);
A = [ones(size(x)), log10(x)];
sol = inv(A'*A)*A'*log10(f);
rate = sol(2)
