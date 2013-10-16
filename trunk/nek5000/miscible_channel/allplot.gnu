! grep ctot channel.log > c.dat

set term postscript enhanced

set xlabel font "Times-Roman, 25"
set ylabel font "Times-Roman, 25"
set xtics font "Times-Roman, 25"
set ytics font "Times-Roman, 25"

set out 'ctot.eps'
set xlabel 'Time'
set ylabel 'M_{0.95}/M_0'
unset key
p 'c.dat' u 5:6 w l lw 2,1-x/40

set out

set out 'xtip.eps'
set xlabel 'Time'
set ylabel 'x_{tip}'
unset key
p 'c.dat' u 5:7 w l lw 2,x
