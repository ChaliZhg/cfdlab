#!/bin/bash

grep Totalke box.log > ke.dat
grep diss_rate box.log > diss.dat
gnuplot plots.gnu
echo "See ke.eps and diss.eps"
