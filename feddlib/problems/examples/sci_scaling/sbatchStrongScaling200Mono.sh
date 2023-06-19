#! /bin/bash -l

#SBATCH -N 8
#SBATCH --ntasks=200
#SBATCH -t 06:30:00
#SBATCH --output=200_Mono.out
#SBATCH --error=200_Mono.err
#SBATCH --switches=1
#SBATCH --cpu-freq=2400000-2400000

unset SLURM_EXPORT_ENV

srun ./problems_sci_scaling.exe

