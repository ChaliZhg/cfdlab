#!/bin/bash

DVPREPRO=$UQTOOL_HOME/utils/dvprepro
PRIMAL=$UQTOOL_HOME/examples/conv_diff/src/flo
ADJOINT=$UQTOOL_HOME/examples/conv_diff/src/adj
CCRE=$UQTOOL_HOME/examples/conv_diff/src/ccre
RECONSTRUCT=$UQTOOL_HOME/examples/conv_diff/src/reconstruct

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

if [ $mode -eq 2 ];
then
   $RECONSTRUCT > recon.log
fi

$PRIMAL > flo.log

$ADJOINT > adj.log

if [ $mode -eq 2 ];
then
   $CCRE > ccre.log
fi
