#!/bin/sh
#PBS -N "pbs.sh"
#PBS -o "mds.log"
#PBS -j oe
#PBS -l "nodes=10:nef:ppn=2"
#PBS -l "walltime=48:00:00"
#PBS -m e
#PBS -q "idle"

cd $PBS_O_WORKDIR
# path to the binary
lamboot
mpirun C ~/work/trunk/optim/src/optim_mpi
lamhalt
