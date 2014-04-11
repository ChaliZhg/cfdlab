set term postscript enhanced
set out 'hist.ps'

set xlabel 't'
set ylabel 'Kinetic energy'
p 'history.dat' u 1:2 t 'Steady KE' w l lw 2, \
  'history.dat' u 1:3 t 'Unsteady KE' w l lw 2

set title 'Evolution of perturbation energy'
set xlabel 't'
set ylabel 'Kinetic energy'
#set logscale y
p 'history.dat' u 1:4 w l lw 2
