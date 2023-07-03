#!/bin/bash -l

#$ -S /bin/bash
#$ -N Index
#$ -j y
#$ -cwd
#$ -pe mpi 1
#$ -q development.q

module load samtools/1.9

samtools index $1
