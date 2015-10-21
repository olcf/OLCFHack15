#PBS -A TRN001
#PBS -l nodes=1,walltime=00:10:00

cd /lustre/atlas/scratch/csep28/trn001/XGC/small_circ
module load cudatoolkit
export PMI_NO_FORK=1

aprun -b -n 1 -N 1 nvprof --events global_ld_mem_divergence_replays,global_st_mem_divergence_replays,gld_request,gst_request /lustre/atlas/scratch/csep28/trn001/XGC/xgc2
export PGI_ACC_TIME=1
aprun -n 1 -N 1 /lustre/atlas/scratch/csep28/trn001/XGC/xgc2

