#!/bin/sh

#SBATCH --time=4:00:00
#SBATCH --job-name=readlog
#SBATCH --partition=cpu_gce
#SBATCH --ntasks=1
#SBATCH --qos=regular
#SBATCH --error="slurm_readlog.err"

#This script is adapted from the example from https://rcc.uchicago.edu/docs/tutorials/kicp-tutorials/running-jobs.html.
#module load parallel

# the --exclusive to srun makes srun use distinct CPUs for each job step
# -N1 -n1 allocates a single core to each task
srun="srun --exclusive -N1 -n1"

# --delay .2 prevents overloading the controlling node
# -j is the number of tasks parallel runs so we set it to $SLURM_NTASKS
# --joblog makes parallel create a log of tasks that it has already run
# --resume makes parallel use the joblog to resume from where it has left off
# the combination of --joblog and --resume allow jobs to be resubmitted if
# necessary and continue from where they left off
#parallel="parallel --delay .2 -j $SLURM_NTASKS --joblog ran_scan_cpv.log --resume"

# this runs the parallel command we want
# in this case, we are running a script named runtask
# parallel uses ::: to separate options. Here {0..99} is a shell expansion
# so parallel will run the command passing the numbers 0 through 99
# via argument {1}
NCPU=30
FILE="outputs/full_potential_cpv/ks210_vs100/scan_r10-4_pt/full_potential/pt"
#$parallel "$srun python ran_scan_cpv_nompi.py {1} $SLURM_NTASKS $FILE $NPT" ::: {0..99}
$srun python test/read_log_pt.py $NCPU $FILE
