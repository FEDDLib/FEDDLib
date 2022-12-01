#! /bin/bash

#SBATCH -J nl
#SBATCH -N 10
#SBATCH --cpus-per-task=24
#SBATCH -t 00:60:00
#SBATCH --output=FSI_240_FaCSITeko_basis.out
#SBATCH --error=FSI_240_FaCSITeko_basis.err

cd $SLURM_SUBMIT_DIR
export I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=0
export nranks=240

#mpirun -np $nranks ./problems_fsi.exe --problemfile=parametersProblemFSI_MonoR_basis.xml --precfileFluidMono=parametersPrecFluidMono_R_basis.xml
#mpirun -np $nranks ./problems_fsi.exe --problemfile=parametersProblemFSI_MonoR.xml --precfileFluidMono=parametersPrecFluidMono_R_full.xml

mpirun -np $nranks ./problems_fsiPaper.exe --problemfile=parametersProblemFSITeko_SIMPLE_R_CB.xml --precfileFluidTeko=parametersPrecFluidTeko_SIMPLE_R_CB.xml
mpirun -np $nranks ./problems_fsiPaper.exe --problemfile=parametersProblemFSITeko_SIMPLE_R.xml --precfileFluidTeko=parametersPrecFluidTeko_SIMPLE_R.xml
