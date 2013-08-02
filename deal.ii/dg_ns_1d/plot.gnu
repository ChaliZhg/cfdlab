set xran[-0.15:-0.11]
p 'solution.gpl' u 1:3 w l, \
  './ns_shock_structure/u' u ($1-0.05645):2 w l

p 'ns.dat' u 1:2 w lp pt 6, \
  './ns_shock_structure/v' u ($1-0.05645):2 w l
