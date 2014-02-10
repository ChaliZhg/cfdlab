filename(n) = sprintf("sol_%d.gpl",n)
N=system("ls -1 sol_*.gpl | wc -l")

set term pdf

set out 'lax_densol.pdf'
set ylabel 'Density'
set xran[-5:5]
#set yran[0.1:1.1]
set key right top
p filename(N-1) u 1:2 t 'DG' w l lw 4, \
  'lax.dat' u ($1-5):2 t 'Exact' w l lw 4 lc 3

set out 'lax_den.pdf'
set ylabel 'Density'
#set yran[0.1:1.1]
set key right top
p 'avg.gpl' u 1:2 t 'DG' w p pt 6 lw 3, \
  'lax.dat' u ($1-5):2 t 'Exact' w l lw 4 lc 3

set out 'lax_vel.pdf'
set ylabel 'Velocity'
#set yran[-0.1:1.0]
set key bottom center
p 'avg.gpl' u 1:3 t 'DG' w p pt 6 lw 3, \
  'lax.dat' u ($1-5):3 t 'Exact' w l lw 4 lc 3

set out 'lax_pre.pdf'
set ylabel 'Pressure'
#set yran[-0.0:1.1]
set key right top
p 'avg.gpl' u 1:4 t 'DG' w p pt 6 lw 3, \
  'lax.dat' u ($1-5):4 t 'Exact' w l lw 4 lc 3

set out 'lax_ind.pdf'
set ylabel 'Indicator'
set yran[-0.0:1.4]
unset key
p 'avg.gpl' u 1:5 t 'DG' w p pt 6 lw 3, \
  'lax.dat' u ($1-5):2 t 'Exact' w l lw 4 lc 3
