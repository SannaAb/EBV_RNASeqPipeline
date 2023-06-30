#!/bin/bash -l

#$ -S /bin/bash
#$ -N TrimGalore
#$ -j y
#$ -cwd
#$ -pe mpi 2
#$ -q development.q

# Running Trimgalore, for removing illumina universal Adapter

echo "Running Trimgalore"

echo "module load trim_galore/0.4.0"
echo "module load cutadapt/1.9"
echo "module load fastqc/0.11.2"

module load trim_galore/0.4.0
module load cutadapt/1.9
module load fastqc/0.11.2

mkdir -p ../Filt_Adaptertrimmed


echo "
trim_galore --fastqc --paired -q 20 --length 30 --output_dir Filt_Adaptertrimmed $1 $2
"

trim_galore --fastqc --paired -q 20 --length 30 --output_dir ../Filt_Adaptertrimmed $1 $2


echo "finished Adaptertriming"



