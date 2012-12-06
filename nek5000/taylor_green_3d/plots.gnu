set term postscript enhanced color
set out 'ke.eps'
set xlabel 'Time'
set ylabel 'Kinetic energy'
p 'ref_data/spectral_Re1600_512.gdiag' u 1:2 t 'Spectral 512^3' w l lw 4, \
  'ke.dat' u 2:3 t 'NEK5000' w l lw 4

set term x11
set out

set term postscript enhanced color
set out 'diss.eps'
set xlabel 'Time'
set ylabel 'Dissipation rate'
p 'ref_data/spectral_Re1600_512.gdiag' u 1:3 t 'Spectral 512^3' w l lw 4, \
  'diss.dat' u 2:3 t 'NEK5000' w l lw 4

set term x11
set out
