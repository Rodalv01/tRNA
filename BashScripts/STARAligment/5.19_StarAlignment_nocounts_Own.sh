#!/usr/bin/env bash

#SBATCH -J StarAlignment_Own_19(nocounts)
#SBATCH --ntasks=1  
#SBATCH --cpus-per-task=6
#SBATCH --mem=32gb              # Amount of RAM needed for this job:
#SBATCH --time=24:00:00         # The time the job will be running:
##SBATCH --gres=gpu:1           # To use GPUs you have to request them:
#SBATCH --array=1-13
#SBATCH -e /mnt/home/users/gm_001_ucm/lsmc/rodrigo/projects/tRNA/logs/slurm_%A_%a.err   #Standard error
#SBATCH -o /mnt/home/users/gm_001_ucm/lsmc/rodrigo/projects/tRNA/logs/slurm_%A_%a.out   #Standard output
##SBATCH --mail-type=END,FAIL
##SBATCH --mail-user=rodalv01@ucm.es
#SBATCH --constraint=cal

######################################################################################################
#Script Name: 5.19_StarAlignment_nocounts_Own.sh 													 #
#Description: Alignment fastq files with reference genome using STAR with starindex(19) that comes   #
#from the tRNA gtf that I create																	 #
#Author: Rodrigo Álvarez Pardo                                          							 #							
#Date: 18/10/2021        																			 #                             
######################################################################################################
set -eu
echo -e "STARTING JOB ${SLURM_ARRAY_TASK_ID}\n"

INDIR="/mnt/home/users/gm_001_ucm/lsmc/rodrigo/projects/tRNA/data/samples/" #Fastq files directory
OUTDIR="/mnt/scratch/users/gm_001_ucm/lsmc/rodrigo/projects/tRNA/output/Star_Alignment/tRNAgtf/3_tRNA.19/2_Anocounts/"  #STAR output directory

SAMPLE_LIST=$(<list_fastq.txt)
SAMPLE_ARRAY=($SAMPLE_LIST)  #Los paréntesis lo convierten en array
SAMPLE=${SAMPLE_ARRAY[${SLURM_ARRAY_TASK_ID}-1]}
echo ${SAMPLE}


genomeOutput="/mnt/home/users/gm_001_ucm/lsmc/rodrigo/projects/tRNA/data/genomeref/STAR_index/GTFtRNA/3_tRNA.19/"


module load star/2.7.9a

STAR --runMode alignReads \
--runThreadN 10 \
--genomeDir ${genomeOutput} \
--readFilesIn /mnt/home/users/gm_001_ucm/lsmc/rodrigo/projects/tRNA/data/samples/${SAMPLE} \
--outFileNamePrefix ${OUTDIR}${SAMPLE} \
--readFilesCommand zcat \
--outStd Log \
--outSAMtype BAM SortedByCoordinate \
--outSAMstrandField intronMotif \
--outFilterMultimapNmax 50 \
--outFilterScoreMinOverLread 0 \
--outFilterMatchNmin 18 \
--outFilterMatchNminOverLread 0 \
--outFilterMismatchNoverLmax 0.04 \
--alignIntronMax 1 

echo -e "\nFINISHED JOB"