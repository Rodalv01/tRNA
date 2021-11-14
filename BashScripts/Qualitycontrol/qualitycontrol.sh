#!/usr/bin/env bash

#SBATCH -J Qualitycontrol
#SBATCH --ntasks=1              # Number of desired cpu:
##SBATCH --cpus-per-task=16
#SBATCH --mem=2gb               # Amount of RAM needed for this job:
#SBATCH --time=24:00:00         # The time the job will be running:
##SBATCH --gres=gpu:1           # To use GPUs you have to request them:
#SBATCH --array=1-13
#SBATCH -e /mnt/home/users/gm_001_ucm/lsmc/rodrigo/projects/tRNA/logs/slurm_%A_%a.err   #Standard error
#SBATCH -o /mnt/home/users/gm_001_ucm/lsmc/rodrigo/projects/tRNA/logs/slurm_%A_%a.out   #Standard output
##SBATCH --mail-type=END,FAIL
##SBATCH --mail-user=rodalv01@ucm.es

############################################################################################################
# Script name: qualitycontrol.sh                                                                           #
# Description: Script to run a quality control test of fasta files with fastqc                             #
# Author: Rodrigo Álvarez Pardo                                                                            #
# Date: 06/09/21                                                                                           #
############################################################################################################

set -eu
echo -e "\nSTARTING JOB\n"

## VARIABLES ###

OUTDIR="/mnt/home/users/gm_001_ucm/lsmc/rodrigo/projects/tRNA/results/"
SAMPLE_LIST=$(<list_fastq.txt)
SAMPLE_ARRAY=($SAMPLE_LIST[@])  #Los paréntesis lo convierten en array
SAMPLE=${SAMPLE_ARRAY[${SLURM_ARRAY_TASK_ID}-1]}
echo ${SAMPLE}

###-----------###

module load fastqc/0.11.4

fastqc -o ${OUTDIR} /mnt/home/users/gm_001_ucm/lsmc/rodrigo/projects/tRNA/data/samples/${SAMPLE} -t 16

echo "FINISHED JOB ${SLURM_ARRAY_TASK_ID}"