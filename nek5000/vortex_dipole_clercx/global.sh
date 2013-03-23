#!/bin/bash

grep Totalke dipole.log > ke.dat
grep enstrophy dipole.log > enstrophy.dat
gnuplot plots.gnu
echo "See ke.eps and enstrophy.eps"
