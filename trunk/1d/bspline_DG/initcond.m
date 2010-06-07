% initial condition
function f = initcond(x)

globals;

n = length(x);

if testcase==lincon_sine
   for j=1:n
      f(j) = sin(2*pi*x(j));
   end
elseif testcase==lincon_hat
   for j=1:n
      if x(j) < 0
         nn = ceil( -x(j) ); x(j) = x(j) + nn;
      end
      if x(j) > 1
         nn = floor(1-x(j)); x(j) = x(j) + nn;
      end
      if x(j) < 0.25 || x(j) > 0.75
         f(j) = 0.0;
      elseif x(j) >= 0.25 && x(j) < 0.5
         f(j) = 4*(x(j) - 0.25);
      else
         f(j) = 1.0 - 4*(x(j) - 0.5);
      end
   end
elseif testcase==burger_sine
   for j=1:n
      f(j) = 0.25 + 0.5 * sin(pi*x(j));
   end
else
   fprintf(1,'Unknown initial condition\n');
   pause
end
