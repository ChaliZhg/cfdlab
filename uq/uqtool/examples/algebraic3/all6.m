Jexact = load('Jexact.dat')

%-----------
Jdata= load('RESULT_P1_uni/J.dat');

Ns= Jdata(:,1);
J = Jdata(:,3);
Jc= Jdata(:,4);

J_error = 100 * abs(J - Jexact) / abs(Jexact);
Jc_error= 100 * abs(Jc - Jexact) / abs(Jexact);

loglog(Ns,J_error,'--r',...
       Ns,Jc_error,'-r',...
       'LineWidth',1.5);

hold on

%richardson extrapolation
for j=2:length(Ns)
   N_rich(j-1) = Ns(j);
   J_rich = (4.0 * J(j) - J(j-1) ) / 3.0;
   err_J_rich(j-1) = 100*abs( (J_rich - Jexact) / Jexact );
end

loglog(N_rich,err_J_rich,'-.r','LineWidth',1.5)

%-----------
Jdata= load('RESULT_P1_adap/J.dat');

Ns= Jdata(:,1);
J = Jdata(:,3);
Jc= Jdata(:,4);

J_error = 100 * abs(J - Jexact) / abs(Jexact);
Jc_error= 100 * abs(Jc - Jexact) / abs(Jexact);

loglog(Ns,J_error,'--b',...
       Ns,Jc_error,'-b',...
       'LineWidth',1.5);

%-----------
Jdata= load('RESULT_P2_uni/J.dat');

Ns= Jdata(:,1);
J = Jdata(:,3);
Jc= Jdata(:,4);

J_error = 100 * abs(J - Jexact) / abs(Jexact);
Jc_error= 100 * abs(Jc - Jexact) / abs(Jexact);

loglog(Ns,J_error,'--m',...
       Ns,Jc_error,'-m',...
       'LineWidth',1.5);

%richardson extrapolation
for j=2:length(Ns)
   N_rich(j-1) = Ns(j);
   J_rich = (8.0 * J(j) - J(j-1) ) / 7.0;
   err_J_rich(j-1) = 100*abs( (J_rich - Jexact) / Jexact );
end

loglog(N_rich,err_J_rich,'-.m','LineWidth',1.5)
%-----------
Jdata= load('RESULT_P2_adap/J.dat');

Ns= Jdata(:,1);
J = Jdata(:,3);
Jc= Jdata(:,4);

J_error = 100 * abs(J - Jexact) / abs(Jexact);
Jc_error= 100 * abs(Jc - Jexact) / abs(Jexact);

loglog(Ns,J_error,'--k',...
       Ns,Jc_error,'-k',...
       'LineWidth',1.5);

set(gca,'FontSize',14)
axis([1 100 1e-7 10])
xlabel('Number of samples')
ylabel('% Error')
legend('P1,uni,J','P1,uni,J+CC','P1,Richardson','P1,adap,J','P1,adap,J+CC','P2,uni,J','P2,uni,J+CC','P2,Richardson','P2,adap,J','P2,adap,J+CC')

hold off
