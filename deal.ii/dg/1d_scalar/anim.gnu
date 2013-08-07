reset
unset key
set grid
set yran[-1.1:1.1]
filename(n) = sprintf("sol_%d.gpl",n)
do for [i=1:100] {
   plot filename(i) u 1:2 w l lw 2
   pause 0.5
}
