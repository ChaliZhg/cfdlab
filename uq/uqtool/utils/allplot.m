% Program to find and plot convergence rate of CC and RE
% Ignore first iskip data points
function allplot ()

Jexact = load('Jexact.dat');

Jdata= load('RESULT/J.dat');
ccre = load('RESULT/error.dat');

J = Jdata(:,3);
Jc= Jdata(:,4);
h = ccre(:,2);
cc= ccre(:,3);
re= ccre(:,4);

J_error = abs(J - Jexact);
Jc_error= abs(Jc - Jexact);

loglog(h,J_error,'o--',...
       h,Jc_error,'*--',...
       h,cc,...
       h,re,'LineWidth',1.5);
legend('Error in J','Error in J+CC','CC','RE')
