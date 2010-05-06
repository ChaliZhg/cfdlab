% initial condition
function f = initcond(x)

globals;

n = length(x);

if testcase==lincon_sine
   for j=1:n
      f(j) = sin(2*pi*x(j));
   end
elseif testcase==burger_sine
   for j=1:n
      f(j) = 0.25 + 0.5 * sin(pi*x(j));
   end
else
   fprintf(1,'Unknown initial condition\n');
   pause
end
