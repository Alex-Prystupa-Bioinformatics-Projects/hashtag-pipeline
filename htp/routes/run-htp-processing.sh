#!/bin/bash
#BSUB -q premium
#BSUB -P acc_naiklab
#BSUB -J run_htp_processing
#BSUB -n 5
#BSUB -R "rusage[mem=30000]"
#BSUB -W 1:00
#BSUB -o logs/run_htp_processing_%J.out
#BSUB -e logs/run_htp_processing_%J.err

module purge
module load R/4.2.0

proj_dir=$(pwd)

mkdir -p logs

echo "RES_COLS=$RES_COLS"

# Pass RES_COLS as an argument to Rscript
bash_cmd="Rscript --vanilla ${proj_dir}/htp/scripts/HTP-Processing.R --res_cols $RES_COLS"

echo $bash_cmd
($bash_cmd)
