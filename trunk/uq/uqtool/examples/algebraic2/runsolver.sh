#!/bin/bash

DVPREPRO=$UQTOOL_HOME/utils/dvprepro
PRIMAL=$UQTOOL_HOME/examples/algebraic2/primal.py
ADJOINT=$UQTOOL_HOME/examples/algebraic2/adjoint.py
RE=$UQTOOL_HOME/examples/algebraic2/re.py

mode=$1
dir=$2
tpl=$3

cp -r $tpl/* $dir
cd $dir

$DVPREPRO random.dat solver.in

$PRIMAL $mode

$ADJOINT $mode

if [ $mode -eq 2 ];
then
   $RE
fi
