set term postscript enhanced color
set out 'ke.eps'
set xlabel 'Time'
set ylabel 'Kinetic energy'
p 'ke.dat' u 2:3 t 'NEK5000' w l lw 4

set term x11
set out

set term postscript enhanced color
set out 'enstrophy.eps'
set xlabel 'Time'
set ylabel 'Enstrophy'
p 'enstrophy.dat' u 2:3 t 'NEK5000' w l lw 4

set term x11
set out
