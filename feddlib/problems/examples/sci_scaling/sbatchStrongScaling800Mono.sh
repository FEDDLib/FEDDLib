#! /bin/bash -l

#SBATCH -N 32
#SBATCH --ntasks=800
#SBATCH -t 06:30:00
#SBATCH --output=800_Mono.out
#SBATCH --error=800_Mono.err
#SBATCH --switches=1
#SBATCH --cpu-freq=2400000-2400000

unset SLURM_EXPORT_ENV

srun ./problems_sci_scaling.exe

