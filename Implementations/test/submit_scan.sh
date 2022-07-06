#!/bin/sh
#SBATCH --job-name=scanbm
#SBATCH --ntasks=100
#SBATCH --partition=cpu_gce
#SBATCH --qos=regular

# This script is adapted from the example script https://rcc.uchicago.edu/docs/tutorials/kicp-tutorials/running-jobs.html.
# Load the default version of GNU parallel.
module purge
module load mpi4py

FILE="outputs/full_potential_cpv/ks23000_vs100/scan_r2"
NPT=1000

mpirun -n $SLURMN_NTASKS python ran_scan_cpv.py $FILE $NPT
