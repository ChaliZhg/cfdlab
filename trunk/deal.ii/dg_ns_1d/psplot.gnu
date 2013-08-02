set term postscript enhanced
set out 'plot.eps'

set xran[-0.15:-0.11]

set ylabel 'Density'
p 'solution.gpl' u 1:2 t 'KFVS' w l, \
  './exact/r' u ($1-0.05645):2 t 'Exact' w l

set ylabel 'Velocity'
p 'solution.gpl' u 1:3 t 'KFVS' w l, \
  './exact/u' u ($1-0.05645):2 t 'Exact' w l

set ylabel 'Pressure'
p 'solution.gpl' u 1:4 t 'KFVS' w l, \
  './exact/p' u ($1-0.05645):2 t 'Exact' w l

set ylabel 'Shear stress'
p 'ns.dat' u 1:2 t 'KFVS' w lp pt 6, \
  './exact/v' u ($1-0.05645):2 t 'Exact' w l

set ylabel 'Heat flux'
p 'ns.dat' u 1:3 t 'KFVS' w lp pt 6, \
  './exact/h' u ($1-0.05645):2 t 'Exact' w l
