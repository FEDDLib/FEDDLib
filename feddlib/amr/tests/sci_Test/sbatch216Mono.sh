#! /bin/bash -l

#SBATCH -N 3
#SBATCH --ntasks=216
#SBATCH -t 03:30:00
#SBATCH --output=216_Mono.out
#SBATCH --error=216_Mono.err
#SBATCH --switches=1
#SBATCH --cpu-freq=2400000-2400000

unset SLURM_EXPORT_ENV

srun ./problems_sci_Test.exe --precfile=parametersPrec_GDSW_CB_CM.xml	--problemfile=parametersProblemSCI_GDSW_CB_CM.xml
srun ./problems_sci_Test.exe --precfile=parametersPrec_RGDSW_CB_CM.xml --problemfile=parametersProblemSCI_RGDSW_CB_CM.xml


