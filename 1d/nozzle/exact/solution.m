function [r u p m] = solution(xs, r1, s1, x)

global gamma frate

a  = nozarea(xs);

% left state
rl = density(xs, s1, r1);
pl = s1*rl^gamma;
ul = frate/(rl*a);
cl = sqrt(gamma*pl/rl);
ml = ul/cl;

if ml>1.0
   % right state
   pr = pl*(1 + 2*gamma*(ml^2-1)/(gamma+1));
   rr = rl*(gamma+1)*ml^2/((gamma-1)*ml^2 + 2);
   ur = ul*rl/rr;
   cr = sqrt(gamma*pr/rr);
   mr = ur/cr;
else
   pr = pl;
   rr = rl;
   ur = ul;
   cr = cl;
   mr = ml;
end

% post-shock entropy
s2 = pr/rr^gamma;

for j=1:length(x)
   a = nozarea(x(j));
   if x(j) < xs
      r(j) = density(x(j), s1, r1);
      p(j) = s1*r(j)^gamma;
   else
      r(j) = density(x(j), s2, rr);
      p(j) = s2*r(j)^gamma;
   end
   u(j) = frate/(r(j)*a);
   c    = sqrt(gamma*p(j)/r(j));
   m(j) = u(j)/c;
end
