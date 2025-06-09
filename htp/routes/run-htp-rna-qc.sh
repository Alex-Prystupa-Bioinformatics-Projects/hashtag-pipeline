#!/bin/bash
#BSUB -q premium
#BSUB -P acc_naiklab
#BSUB -J run_htp
#BSUB -n 5
#BSUB -R "rusage[mem=30000]"
#BSUB -W 1:00
#BSUB -o logs/run_rna_qc_%J.out
#BSUB -e logs/run_rna_qc_%J.err

module purge
module load R/4.2.0

proj_dir=$(pwd)

mkdir logs

# bash commands
bash_cmd="Rscript --vanilla ${proj_dir}/htp/scripts/RNA-QC.R"

# Run Step-1
echo $bash_cmd
($bash_cmd) 
