#! /bin/bash -l

#SBATCH -N 16
#SBATCH --ntasks=400
#SBATCH -t 06:30:00
#SBATCH --output=400_Mono.out
#SBATCH --error=400_Mono.err
#SBATCH --switches=1
#SBATCH --cpu-freq=2400000-2400000

unset SLURM_EXPORT_ENV

srun ./problems_sci_scaling.exe

