#! /bin/bash -l

#SBATCH -N 4
#SBATCH --ntasks=100
#SBATCH -t 06:30:00
#SBATCH --output=100_Mono.out
#SBATCH --error=100_Mono.err
#SBATCH --switches=1
#SBATCH --cpu-freq=2400000-2400000

unset SLURM_EXPORT_ENV

srun ./problems_sci_scaling.exe

