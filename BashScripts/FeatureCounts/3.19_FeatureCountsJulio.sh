#!/usr/bin/env bash

#SBATCH -J FeatureCount_Julio(v19)
#SBATCH --ntasks=1  
#SBATCH --cpus-per-task=6
#SBATCH --mem=32gb              # Amount of RAM needed for this job:
#SBATCH --time=24:00:00         # The time the job will be running:
##SBATCH --gres=gpu:1           # To use GPUs you have to request them:
#SBATCH --array=1
#SBATCH -e /mnt/home/users/gm_001_ucm/lsmc/rodrigo/projects/tRNA/logs/slurm_%A_%a.err   #Standard error
#SBATCH -o /mnt/home/users/gm_001_ucm/lsmc/rodrigo/projects/tRNA/logs/slurm_%A_%a.out   #Standard output
##SBATCH --mail-type=END,FAIL
##SBATCH --mail-user=rodalv01@ucm.es
#SBATCH --constraint=cal

######################################################################################################
#Script Name: 3.19_FeatureCountsJulio.sh 														 	 #
#Description: Alignment fastq files with reference genome using STAR with starindex that comes from  #
# Julio gtf																						     #
#Author: Rodrigo Álvarez Pardo                                          							 #							
#Date: 13/10/2021        																			 #                             
######################################################################################################

GTFDIR="/mnt/home/users/gm_001_ucm/lsmc/rodrigo/projects/tRNA/data/annotation/GTFs_tRNAs/"
BAMDIR="/mnt/scratch/users/gm_001_ucm/lsmc/rodrigo/projects/tRNA/output/Star_Alignment/Juliogtf/3_Juliogtf.19/2_Anocounts/*.bam"
OUTDIR="/mnt/home/users/gm_001_ucm/lsmc/rodrigo/projects/tRNA/outputs/Feature_count/1_Juliogtf/3_Featurecounts.19/"

module load subread

featureCounts -a ${GTFDIR}Spiked_GRCm38_tRNAs.gtf \
-o ${OUTDIR}Featurecounts.count ${BAMDIR} -M --fraction -s 2 -T 6

echo "Finish Job"