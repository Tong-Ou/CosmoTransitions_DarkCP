#!/bin/sh

#SBATCH --time=24:00:00
#SBATCH --job-name=calcpt
#SBATCH --partition=cpu_gce
#SBATCH --ntasks=160
#SBATCH --qos=regular

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
OUTPATH="outputs/full_potential_cpv/wrap_up/ma10to350_rs1e3_pt/full_potential/"
mkdir -p $OUTPATH
#$parallel "$srun python pt_cpv.py {1} $SLURM_NTASKS $FILE $OUTPATH" ::: {0..99}
for i in {0..159}
do
        $srun python pt_forGW.py $i $SLURM_NTASKS $OUTPATH&
done
wait
