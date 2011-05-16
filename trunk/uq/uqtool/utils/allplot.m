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

J_error = 100 * abs(J - Jexact) / abs(Jexact);
Jc_error= 100 * abs(Jc - Jexact) / abs(Jexact);
cc = 100 * cc / abs(Jexact);
re = 100 * re / abs(Jexact);

figure(1)
plot(Ns,J,'o-',Ns,Jc,'*--',Ns,ones(size(Ns))*Jexact,'--','LineWidth',1.5)
set(gca,'FontSize',14)
xlabel('Number of samples')
ylabel('Mean functional')
legend('J','J+CC','Exact')

figure(2)
loglog(Ns,J_error,'o--',...
       Ns,Jc_error,'*--',...
       Ns,cc,...
       Ns,re,'LineWidth',1.5);
set(gca,'FontSize',14)
xlabel('Number of samples')
ylabel('Percentage error')
legend('Error in J','Error in J+CC','CC','RE')

% Convergence rates
ns = length(h);
x  = Ns(istart:ns);

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
