#! /bin/bash -l

#SBATCH -N 64
#SBATCH --ntasks=4096
#SBATCH -t 04:30:00
#SBATCH --output=4096_Mono.out
#SBATCH --error=4096_Mono.err
#SBATCH --switches=1
#SBATCH --cpu-freq=2400000-2400000

unset SLURM_EXPORT_ENV

srun ./problems_sci_scaling.exe 

