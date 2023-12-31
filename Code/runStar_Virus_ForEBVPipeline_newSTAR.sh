#!/bin/bash -l

#$ -N VirusAlign_EBV
#$ -j y
#$ -cwd
#$ -pe mpi 40
#$ -q development.q

#genome=/medstore/projects/B22-017/Intermediate/db/star_2.7.2b_NC_007605.1_index
genome=""

#STAR=~/Programs/STAR-2.7.2b/bin/Linux_x86_64/STAR
STAR=""

mkdir -p ../EBV_Alignment

outname=../EBV_Alignment/${1%_R1_001_val_1.fq.gz}_NC_007605.1


echo "
$STAR --runThreadN 40 --genomeDir $genome --readFilesIn $1 $2 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $outname --limitBAMsortRAM 100000000000 --chimOutType WithinBAM --outSAMunmapped Within --outFilterMultimapNmax 100 --readFilesCommand zcat"


$STAR --runThreadN 40 --genomeDir $genome --readFilesIn $1 $2 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $outname --limitBAMsortRAM 100000000000 --chimOutType WithinBAM --outSAMunmapped Within --outFilterMultimapNmax 100 --readFilesCommand zcat


echo "Fin"
