Jexact = load('Jexact.dat')

Jdata= load('RESULT_eno_adap/J.dat');

Ns= Jdata(:,1);
J = Jdata(:,3);
Jc= Jdata(:,4);

J_error = 100 * abs(J - Jexact) / abs(Jexact);
Jc_error= 100 * abs(Jc - Jexact) / abs(Jexact);

figure(1)
plot(Ns,J,'o-',Ns,Jc,'*--',Ns,ones(size(Ns))*Jexact,'-.','LineWidth',1.5)
set(gca,'FontSize',14)
axis([0 70 3 3.3])
xlabel('Number of samples')
ylabel('Mean functional')
legend('J','J+CC','Exact')

figure(2)
loglog(Ns,J_error,'o-',...
       Ns,Jc_error,'*-',...
       'LineWidth',1.5);
set(gca,'FontSize',14)
axis([1 100 1e-7 10])
xlabel('Number of samples')
ylabel('% Error')
legend('J','J+CC')
