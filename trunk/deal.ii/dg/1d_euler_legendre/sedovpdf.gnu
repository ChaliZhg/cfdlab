set grid
filename(n) = sprintf("sol_%d.gpl",n)
N=system("ls -1 sol_*.gpl | wc -l")

set term pdf

set out 'sed_densol.pdf'
set key center
set yran[-0.1:6]
set ylabel 'Density'
set xlabel 'x'
set title 'DG Solution: Density'
plot 'sedov.dat' u 2:3 t 'Exact' w l lw 4 lt 2 lc 3, \
     'sedov.dat' u (-$2):3 t '' w l lw 4 lt 2 lc 3, \
     filename(N-1) u 1:2 t 'DG' w l lt 1 lw 4

set out 'sed_den.pdf'
set key center
set yran[-0.1:6]
set ylabel 'Density'
set xlabel 'x'
set title 'Cell average'
plot 'sedov.dat' u 2:3 t 'Exact' w l lw 4 lt 2 lc 3, \
     'sedov.dat' u (-$2):3 t '' w l lw 4 lt 2 lc 3, \
     'avg.gpl' u 1:2 t 'DG' w p pt 6 lc 1

set out 'sed_vel.pdf'
set key left top
set yran[-810:810]
set ylabel 'Velocity'
set xlabel 'x'
set title 'Cell average'
plot 'sedov.dat' u 2:6 t 'Exact' w l lw 4 lt 2 lc 3, \
     'sedov.dat' u (-$2):(-$6) t '' w l lw 4 lt 2 lc 3, \
     'avg.gpl' u 1:3 t 'DG' w p pt 6 lc 1

set out 'sed_pre.pdf'
set key center
set yran[-0.1:800000]
set ylabel 'Pressure' offset -2,0
set xlabel 'x'
set title 'Cell average'
plot 'sedov.dat' u 2:5 t 'Exact' w l lw 4 lt 2 lc 3, \
     'sedov.dat' u (-$2):5 t '' w l lw 4 lt 2 lc 3, \
     'avg.gpl' u 1:4 t 'DG' w p pt 6 lc 1

set out 'sed_ind.pdf'
set ylabel 'Indicator'
set yran[-0.0:6]
p 'avg.gpl' u 1:5 t 'Indicator' w p pt 6 lw 3, \
  'sedov.dat' u 2:3 t 'Density' w l lw 4 lt 2 lc 3, \
  'sedov.dat' u (-$2):3 t '' w l lw 4 lt 2 lc 3
