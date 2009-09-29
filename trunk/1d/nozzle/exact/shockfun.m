function f = shockfun(xs, x2, r1, s1, p2)

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

% outflow conditions
r2 = density(x2, s2, rr);
pc = s2*r2^gamma;

% function
f = abs(pc - p2)/p2;
