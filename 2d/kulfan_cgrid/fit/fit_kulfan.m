clear all

N1 = 0.5;
N2 = 1.0;

degu = 8;
degl = 8;

load airfoil.txt

% Total number of points
% np MUST BE ODD
np = length(airfoil(:,1));

if mod(np,2) ~= 1
   fprintf(1,'Number of airfoil points must be odd\n');
   pause
end

nu = (np+1)/2;
nl = (np+1)/2;

xl = airfoil(nl:-1:1,1);
yl = airfoil(nl:-1:1,2);

xu = airfoil(nu:np,1);
yu = airfoil(nu:np,2);

if min(xu) == min(xl)
   xle = min(xu);
else
   fprintf(1,'LE point different on upper and lower surface\n');
   pause
end

plot(xu,yu,'o-',xl,yl,'*-')

% Check if first point on upper and lower curve are same
if xl(1) ~= xu(1) || yl(1) ~= yu(1)
   fprintf('LE points on lower and upper surface do not match\n');
   pause
end

% If LE does not have y=0, then shift up/down
if yl(1) ~= 0.0
   yle= yl(1);
   yl = yl - yle;
   yu = yu - yle;
end

if yl(nl) ~= yu(nu) % Then TE is not sharp

   % St line approximation at TE
   % Eqn of upper surface: y = mu*(x - xu(nu)) + yu(nu)
   %                       mu= (yu(nu) - yu(nu-1))/(xu(nu) - xu(nu-1))
   % Eqn of lower surface: y = ml*(x - xl(nl)) + yl(nl)
   %                       ml= (yl(nl) - yl(nl-1))/(xl(nl) - xl(nl-1))

   mu= (yu(nu) - yu(nu-1))/(xu(nu) - xu(nu-1));
   ml= (yl(nl) - yl(nl-1))/(xl(nl) - xl(nl-1));

   % X-coordinate of TE found by intersection of two st. lines
   xte = (mu*xu(nu) - ml*xl(nl) + yl(nl) - yu(nu) )/(mu - ml);
   yte = mu*(xte - xu(nu)) + yu(nu);
   fprintf('xte = %f,   yte = %f\n', xte, yte);

   % Append TE point
   xu = [xu; xte];
   yu = [yu; yte];
   xl = [xl; xte];
   yl = [yl; yte];

end

% Update size
nu = length(xu);
nl = length(xl);

% Translate so that xle = 0
xu = xu - xle;
xl = xl - xle;
xte= xte - xle;
xle= 0.0;

% Airfoil chord = max(xl) = max(xu)
chord = max(xl);
fprintf('chord = %f\n', chord)

plot(xu,yu,'o',xl,yl,'*')
axis tight
axis equal

% Fit kulfan representation

psil = xl/chord;
psiu = xu/chord;

zetal = yl/chord;
zetau = yu/chord;

ztl   = zetal(nl);
ztu   = zetau(nu);

% Shape function evaluated at given points
% Remove LE and TE points since shape function  not defined
Sl = (zetal(2:nl-1) - psil(2:nl-1)*ztl) ./ (psil(2:nl-1)).^N1 ./ ...
     (1-psil(2:nl-1)).^N2;
Su = (zetau(2:nu-1) - psiu(2:nu-1)*ztu) ./ (psiu(2:nu-1)).^N1 ./ ...
     (1-psiu(2:nu-1)).^N2;

% Lower curve
Al = zeros(nl-2,degl+1);
for j=2:nl-1
   for k=0:degl
      Al(j-1,k+1) = bernstein(degl,k,psil(j));
   end
end
Cl = inv(Al' * Al) * Al' * Sl;
errl = sqrt(sum((Al*Cl - Sl).^2)/(nl-2));
fprintf('Shape fun lower surface error = %e\n', errl)

% Upper curve
Au = zeros(nu-2,degu+1);
for j=2:nu-1
   for k=0:degu
      Au(j-1,k+1) = bernstein(degu,k,psiu(j));
   end
end
Cu = inv(Au' * Au) * Au' * Su;
erru = sqrt(sum((Au*Cu - Su).^2)/(nu-2));
fprintf('Shape fun upper surface error = %e\n', erru)

% Plot airfoil using kulfan
x = linspace(0,pi,200);
x = 0.5*(1 - cos(x));
for j=1:length(x)
   shape = decas(Cl, x(j));
   ynl(j) = shape * (x(j))^N1 * (1-x(j))^N2 + x(j)*ztl;
   ynl(j) = ynl(j) * chord;
   xnl(j) = x(j) * chord;

   shape = decas(Cu, x(j));
   ynu(j) = shape * (x(j))^N1 * (1-x(j))^N2 + x(j)*ztu;
   ynu(j) = ynu(j) * chord;
   xnu(j) = x(j) * chord;
end

hold on
plot(xnl,ynl,'-',xnu,ynu,'-')

% Save kulfan into file
fid = fopen('kulfan.out','w');
fprintf(fid,'%12.4f    # N1\n', N1);
fprintf(fid,'%12.4f    # N2\n', N2);
fprintf(fid,'%8d       # degree of bezier curve\n', degl);
fprintf(fid,'# read bezier coefficients for lower surface\n');
fprintf(fid,'%24.14e\n', Cl);
fprintf(fid,'%24.14e   # zite_l\n', ztl);
fprintf(fid,'# read bezier coefficients for upper surface\n');
fprintf(fid,'%24.14e\n', Cu);
fprintf(fid,'%24.14e   # zite_u\n', ztu);
fclose(fid);
