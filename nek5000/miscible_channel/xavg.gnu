reset
set term gif animate
set out 'xavg.gif'
set grid
set xlabel 'ux, c'
set ylabel 'y'
set xran[0:2]
set key autotitle columnhead
N=system("ls -1 channel0.f* | wc -l")
do for [n=1:N] {
    plot 'xavg.dat' i n-1 u 2:1  t columnheader(1) w l lw 2, \
         'xavg.dat' i n-1 u 3:1  t columnheader(2) w l lw 2
}

reset
set term gif animate
set out 'muavg.gif'
set grid
unset key
set xlabel 'mu'
set ylabel 'y'
set xran[0.001:0.02]
N=system("ls -1 channel0.f* | wc -l")
do for [n=1:N] {
    plot 'xavg.dat' i n-1 u 4:1  w l lw 2
}

reset
set term gif animate
set out 'tavg.gif'
set grid
unset key
set xlabel 'T'
set ylabel 'y'
#set xran[0.001:0.02]
N=system("ls -1 channel0.f* | wc -l")
do for [n=1:N] {
    plot 'xavg.dat' i n-1 u 5:1  w l lw 2
}
