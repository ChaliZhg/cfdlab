#!/bin/bash

DVPREPRO=$UQTOOL_HOME/utils/dvprepro
PRIMAL=$UQTOOL_HOME/examples/conv_diff/src/flo
ADJOINT=$UQTOOL_HOME/examples/conv_diff/src/adj
RE=$UQTOOL_HOME/examples/conv_diff/src/ccre

mode=$1
dir=$2
tpl=$3

cp -r $tpl/* $dir
cd $dir

if [ $mode -eq 1 ];
then
   cp param.in.1 param.in
else
   cp param.in.2 param.in
fi

$DVPREPRO random.dat param.in

$PRIMAL > flo.log

$ADJOINT > adj.log

if [ $mode -eq 2 ];
then
   $RE
fi
