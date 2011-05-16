#!/bin/bash

tpl=$1

REFINE=$UQTOOL_HOME/examples/conv_diff/src/refine

cd $tpl
$REFINE > refine.log
