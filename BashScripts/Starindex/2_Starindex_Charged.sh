#!/usr/bin/env bash

#SBATCH -J Starindex_Charged
#SBATCH --ntasks=1              # Number of desired cpu:
#SBATCH --cpus-per-task=6
#SBATCH --mem=32gb              # Amount of RAM needed for this job:
#SBATCH --time=24:00:00         # The time the job will be running:
##SBATCH --gres=gpu:1           # To use GPUs you have to request them:
#SBATCH --array=1
#SBATCH -e /mnt/home/users/gm_001_ucm/lsmc/rodrigo/projects/tRNA/logs/slurm_%A_%a.err   #Standard error
#SBATCH -o /mnt/home/users/gm_001_ucm/lsmc/rodrigo/projects/tRNA/logs/slurm_%A_%a.out   #Standard output
##SBATCH --mail-type=END,FAIL
##SBATCH --mail-user=rodalv01@ucm.es

###################################################################################################
# Script Name: Starindex_Charged.sh                                                                      #
# Description: Script to create the alignment index with start, it is necesary a reference genome #
# (fasta file) which is GRCm38_tRNA.p6 with the tRNA fasta file at the end and a GTF with the tRNA#                                         
# added by me 																					  #
# Author: Rodrigo Álvarez Pardo																	  #
# Date: 18/10/21																				  #
###################################################################################################

set -eu
echo -e "\nSTARTING JOB\n"

#### VARIABLES ####

INDIRG="/mnt/home/users/gm_001_ucm/lsmc/rodrigo/projects/tRNA/data/genomeref/"
echo $INDIRG
INDIRGT="/mnt/home/users/gm_001_ucm/lsmc/rodrigo/projects/tRNA/data/annotation/GTFs_tRNAs/"
echo $INDIRGT

###################

module load star/2.7.9a

STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir /mnt/home/users/gm_001_ucm/lsmc/rodrigo/projects/tRNA/data/genomeref/STAR_index/GTFCharged \
--genomeFastaFiles ${INDIRG}Mus_musculus.GRCm38.dna.primary_assembly.fa \
--sjdbGTFfile ${INDIRGT}Charged_Mouse_tRNAs.gtf

echo -e "\nFINISHED JOB\n"