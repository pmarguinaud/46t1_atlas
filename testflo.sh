#!/bin/bash
#SBATCH -N1
#SBATCH --time 00:05:00
#SBATCH --export="NONE"

set -x

ulimit -s unlimited
ulimit -l unlimited

unset DR_HOOK_NOT_MPI
export DR_HOOK=1
export DR_HOOK_OPT=prof
export DR_HOOK_SILENT=1
export EC_MEMINFO=0
export DR_HOOK_SHOW_PROCESS_OPTIONS=0
export DR_HOOK_IGNORE_SIGNALS=-1
export MPL_MBX_SIZE=2000000000
export OMP_STACKSIZE=4G
NN=$SLURM_NNODES
NNP=8
let "NP=$NN*$NNP"
#export MPIAUTOCONFIG=/home/gmap/mrpa/suzat/.mpiautorc/mpiauto.DDT.conf


cd $TMPDIR

/opt/softs/mpiauto/mpiauto --wrap --wrap-stdeo \
  --prefix-mpirun "/usr/bin/time -f 'real=%e'" \
  --verbose -nn $NN -nnp $NNP -openmp 1 --  /home/gmap/mrpm/marguina/pack/46t1_atlas.01.I185274INTELMPI184274MT.x/bin/TESTFLO

ls -lrt
