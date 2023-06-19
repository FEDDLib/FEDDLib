#! /bin/bash -l

#SBATCH -N 3
#SBATCH --ntasks=216
#SBATCH -t 03:30:00
#SBATCH --output=216_Mono.out
#SBATCH --error=216_Mono.err
#SBATCH --switches=1
#SBATCH --cpu-freq=2400000-2400000

unset SLURM_EXPORT_ENV

srun ./problems_sci_scaling.exe 

