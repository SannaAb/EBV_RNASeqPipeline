#!/bin/bash -l

#$ -N STAR_Human
#$ -j y
#$ -cwd
#$ -pe mpi 40
#$ -q development.q



#genome=/medstore/projects/B22-017/Intermediate/db/StarIndex_2.7.2b_Grch38_index/
genome=""

#STAR=~/Programs/STAR-2.7.2b/bin/Linux_x86_64/STAR
STAR=""

mkdir -p ../Alignment_Grch38

outname=../Alignment_Grch38/${1%_R1_001_val_1.fq.gz}_Grch38

echo "
$STAR --runThreadN 40 --genomeDir $genome --readFilesIn $1 $2 --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --outFileNamePrefix $outname --limitBAMsortRAM 10000000000 --chimOutType WithinBAM 
"

$STAR --runThreadN 40 --genomeDir $genome --readFilesIn $1 $2 --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --outFileNamePrefix $outname --limitBAMsortRAM 10000000000 --chimOutType WithinBAM 

echo "Finished"


