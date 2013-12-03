reset
set grid
filename(n) = sprintf("sol_%d.gpl",n)
N=system("ls -1 sol_*.gpl | wc -l")

set term postscript enhanced color
set out 'result.eps'

set key center
set yran[-0.1:6]
set ylabel 'Density'
set xlabel 'x'
set title 'DG Solution: Density'
plot 'sedov.dat' u 2:3 t 'Exact' w l lw 2 lt 2 lc 3, \
     'sedov.dat' u (-$2):3 t '' w l lw 2 lt 2 lc 3, \
     filename(N-1) u 1:2 t 'DG' w l lt 1 lw 2

set key center
set yran[-0.1:6]
set ylabel 'Density'
set xlabel 'x'
set title 'Cell average'
plot 'sedov.dat' u 2:3 t 'Exact' w l lw 2 lt 2 lc 3, \
     'sedov.dat' u (-$2):3 t '' w l lw 2 lt 2 lc 3, \
     'avg.gpl' u 1:2 t 'DG' w p pt 6 lc 1

set key left top
set yran[-810:810]
set ylabel 'Velocity'
set xlabel 'x'
set title 'Cell average'
plot 'sedov.dat' u 2:6 t 'Exact' w l lw 2 lt 2 lc 3, \
     'sedov.dat' u (-$2):(-$6) t '' w l lw 2 lt 2 lc 3, \
     'avg.gpl' u 1:3 t 'DG' w p pt 6 lc 1

set key center
set yran[-0.1:800000]
set ylabel 'Pressure'
set xlabel 'x'
set title 'Cell average'
plot 'sedov.dat' u 2:5 t 'Exact' w l lw 2 lt 2 lc 3, \
     'sedov.dat' u (-$2):5 t '' w l lw 2 lt 2 lc 3, \
     'avg.gpl' u 1:4 t 'DG' w p pt 6 lc 1
