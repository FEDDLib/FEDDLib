#! /bin/bash -l

#SBATCH -N 14
#SBATCH --ntasks=1000
#SBATCH -t 03:30:00
#SBATCH --output=1000_Mono.out
#SBATCH --error=1000_Mono.err
#SBATCH --switches=1
#SBATCH --cpu-freq=2400000-2400000

unset SLURM_EXPORT_ENV

srun ./problems_sci_scaling.exe

