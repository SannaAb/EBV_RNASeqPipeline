#!/bin/bash -l
#$ -N Calmd
#$ -j y
#$ -cwd
#$ -pe mpi 1
#$ -q development.q

# We need to do this as the new star version does not report the NM which we need in the bamfile for the EBV pipeline

module load samtools/1.9


#genome=/medstore/projects/B22-017/Intermediate/db/fasta/NC_007605.1.fasta

genome=""

echo "

samtools calmd $1 $genome | samtools view -Sb - > ${1%.bam}_AddFlags.bam

"

samtools calmd $1 $genome | samtools view -Sb - > ${1%.bam}_AddFlags.bam

echo "Fin"
