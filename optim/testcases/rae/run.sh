#!/bin/bash

DEFORM=$HOME/work/euler2d/src-mesh/deform
SOLVER=$HOME/work/euler2d/src-flo/flo
MODE=$1
DIR=$2

if [ $MODE == "deform" ]
then
   cd $DIR
   $DEFORM ../flo2.inp > flo.log
fi

if [ $MODE == "solve" ]
then
   cd $DIR
   echo "0 0.0 1.0 1.0e20" > FLO.OUT
   cp ../flo.gnu .
   $SOLVER ../flo.inp >> flo.log
fi
