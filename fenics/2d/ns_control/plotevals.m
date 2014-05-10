load evals50.mat
d1 = D;
load evals100.mat
d2 = D;

plot(real(d1),imag(d1),'o',real(d2),imag(d2),'*','LineWidth',1.5)
legend('50x50','100x100')
grid on
xlabel('real(\lambda)','FontSize',14)
ylabel('imag(\lambda)')
set(gca,'FontSize',16)
