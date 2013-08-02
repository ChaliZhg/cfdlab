function dy=Fb(t,y)
%y(1): T
%y(2): -U
%t: -x 

M=1.5; 
gamma=5/3; 
Pr=2/3; 

%upstream states: 
r1=1.0;
u1=1.0; 
mu1=0.0005; 
p1=1/gamma/M^2; 
T1=2*p1/r1; 
A=r1*u1; 
B=r1*u1^2+p1; 
C=(1/2*r1*u1^2+p1/(gamma-1)+p1)*u1; 
mu=mu1*(y(1)/T1)^0.8; 

%mu=mu1; 
dy=[0;0]; 
dy(1)=4/5*Pr/mu*(1/2*A*y(2)^2-3/4*A*y(1)+C+B*y(2)); % -T x 
dy(2)=3/4/mu*(-B-A*y(2)-1/2*A*y(1)/y(2)); % U x
