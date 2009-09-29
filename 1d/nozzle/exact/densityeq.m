function f = densityeq(r, x, s1)

global gamma frate H0

a = nozarea(x);
f = gamma*s1*r^(gamma-1)/(gamma-1) + 0.5*frate^2/(r^2*a^2) - H0;
