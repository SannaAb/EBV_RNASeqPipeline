#!/bin/bash -l
#$ -N EBVMasterPipe
#$ -j y
#$ -cwd
#$ -pe mpi 1
#$ -q development.q


module load miniconda/4.14.0
#source activate py2
source activate /home/xabras/.conda/envs/py2

# Loop in the human folder!

module load samtools/1.9
module load igvtools/2.1.7

PathToEBVAlignments=""
PathToHumanAlignment=""
PathToFilteredFastq=""

#PathToEBVMasterpipeline=EBV_Master_pipeline3.py
PathToEBVMasterpipeline=""

echo "

python $PathToEBVMasterpipeline -IH $PathToHumanAlignment/$1_Grch38Aligned.sortedByCoord.out_AddFlags.bam -IE $PathToEBVAlignments/$1_NC_007605.1Aligned.sortedByCoord.out_AddFlags.bam -O $1_KeptEBER -F $PathToFilteredFastq/$1_R1_001_val_1.fq.gz -MM 3 -AL 40 -MU 10

"

python $PathToEBVMasterpipeline  -IH $PathToHumanAlignment/$1_Grch38Aligned.sortedByCoord.out_AddFlags.bam -IE $PathToEBVAlignments/$1_NC_007605.1Aligned.sortedByCoord.out_AddFlags.bam -O $1_KeptEBER -F $PathToFilteredFastq/$1_R1_001_val_1.fq.gz -MM 3 -AL 40 -MU 10


echo "Finished"
