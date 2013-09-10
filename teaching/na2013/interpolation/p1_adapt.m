clear all

xmin = 0;
xmax = 1;
fun = @(x) exp(-100*(x-0.5).^2);

N = 10;
x = linspace(xmin,xmax,N);
f = fun(x);

ne = 100;
xe = linspace(xmin,xmax,ne);
fe = fun(xe);
plot(xe,fe,'-',x,f,'o-','LineWidth',2);
title('Initial approximation')
pause

err = 1e-3;
nadapt = 50;

for n=1:nadapt
   h = x(2:end) - x(1:end-1);
   D = zeros(1,N);
   D(2:end-1) = (f(3:end) - f(1:end-2)) ./ (x(3:end) - x(1:end-2));
   H(1:N-1) = (D(2:end) - D(1:end-1)) ./ h;
   elem_err = (1.0/8.0) * h.^2 .* abs(H);
   [current_err,i] = max(elem_err);
   if current_err < err
      fprintf('Satisfied error tolerance\n')
      break
   end
   x = [x(1:i), 0.5*(x(i)+x(i+1)), x(i+1:end)];
   f = [f(1:i), fun(x(i+1)), f(i+1:end)];
   plot(xe,fe,'-',x,f,'or--','LineWidth',2);
   N = length(x);
   pause
end
