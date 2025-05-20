#!/bin/bash
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=Aleksandr.Prystupa@nyulangone.org
#SBATCH --time=01:00:00
#SBATCH --partition=a100_short
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=64G
#SBATCH --job-name=run_htp
#SBATCH --output=logs/run_rna_qc_%A_%a.out
#SBATCH --error=logs/run_rna_qc_%A_%a.err

module purge
module load r/4.2.2

proj_dir=$(pwd)

mkdir logs

# bash commands
bash_cmd="Rscript --vanilla ${proj_dir}/htp/scripts/RNA-QC.R"

# Run Step-1
echo $bash_cmd
($bash_cmd) 
