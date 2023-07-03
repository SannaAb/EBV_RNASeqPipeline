#!/bin/bash -l

#$ -N Calmd
#$ -j y
#$ -cwd
#$ -pe mpi 1
#$ -q development.q

# We need to do this as the new star version does not report the NM which we need in the bamfile for the EBV pipeline

module load samtools/1.9


#genome=/medstore/databases/research/tmp/Homo_sapiens/Ensembl/GRCh38.90/WholeGenomeFasta/Homo_sapiens.GRCh38.dna.toplevel.canonical.fa
genome=""

echo "

samtools calmd $1 $genome | samtools view -Sb - > ${1%.bam}_AddFlags.bam

"

samtools calmd $1 $genome | samtools view -Sb - > ${1%.bam}_AddFlags.bam

echo "Fin"
