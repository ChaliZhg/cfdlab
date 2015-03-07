% y should be in [0.5, 1.0]
function c = c0(y)

h = 0.3;
q = 0.05;

n = length(y);
c = zeros(size(y));

yy = 2*(y - 0.5);

for j=1:n
   if yy(j) < 2.0*h
      c(j) = 0.0;
   elseif yy(j) > 2.0*(h+q)
      c(j) = 1.0;
   else
      xi = (yy(j) - 2*h)/(2*q);
      c(j) = 30 * xi^3 * (1.0/3.0 - 0.5*xi + xi^2/5.0);
   end
end
