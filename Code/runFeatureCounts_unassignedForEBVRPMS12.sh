#!/bin/bash -l

#$ -N FeatureCounts
#$ -j y
#$ -cwd
#$ -pe mpi 1
##$ -q development.q

module load samtools/1.9


#FirstAnno=NC_007605.1-20210315_NC_changedfornon-strandedsequencing.tsv
FirstAnno=""
#SecondAnno=NC_007605.1-20210304_RPMS1unassigned.tsv
SecondAnno=""

#featureCounts=/home/xabras/Programs/subread-2.0.0-source/bin/featureCounts
featureCounts=""

#MergeFeatureCountsscript=/medstore/projects/B22-017/Intermediate/90-843202466/Code/MergeCountsFromTwoFeatureCountsRuns.py
MergeFeatureCountsscript=""


mkdir -p NewCounts




echo "

$featureCounts -a $FirstAnno -M -s 0 -J -p -t exon --Rpath NewCounts/ -R BAM -g gene -o NewCounts/$1.counts EBV_Alignment/$1_KeptEBER.3Missmatch_40nt_10Multimapped.EBV.bam

"

$featureCounts -a $FirstAnno -M -s 0 -J -p -t exon --Rpath NewCounts/ -R BAM -g gene -o NewCounts/$1.counts EBV_Alignment/$1_KeptEBER.3Missmatch_40nt_10Multimapped.EBV.bam


echo "

samtools view NewCounts/$1_KeptEBER.3Missmatch_40nt_10Multimapped.EBV.bam.featureCounts.bam -h | grep -e "^@" -e Unassigned_NoFeatures  | samtools view -Sb - > NewCounts/$1_unassigned.bam

"

samtools view NewCounts/$1_KeptEBER.3Missmatch_40nt_10Multimapped.EBV.bam.featureCounts.bam -h | grep -e "^@" -e Unassigned_NoFeatures  | samtools view -Sb - > NewCounts/$1_unassigned.bam

echo "

samtools sort NewCounts/$1_unassigned.bam -o NewCounts/$1_unassigned_sorted.bam

"

samtools sort NewCounts/$1_unassigned.bam -o NewCounts/$1_unassigned_sorted.bam

echo "

samtools index NewCounts/$1_unassigned_sorted.bam

"

samtools index NewCounts/$1_unassigned_sorted.bam


echo "

$featureCounts -a $SecondAnno -M -s 0 -J -p -t exon -g gene -o NewCounts/$1.Genecounts NewCounts/$1_unassigned_sorted.bam

"

$featureCounts -a $SecondAnno -M -s 0 -J -p -t exon -g gene -o NewCounts/$1.Genecounts NewCounts/$1_unassigned_sorted.bam

echo "

python $MergeFeatureCountsscript NewCounts/$1.counts NewCounts/$1.Genecounts

"


python $MergeFeatureCountsscript NewCounts/$1.counts NewCounts/$1.Genecounts

echo "

cat NewCounts/$1_Merge_unassigned_FeatureCounts.txt Human_counts/$1_count > NewCounts/$1_HV.count

"

#cat NewCounts/$1_Merge_unassigned_FeatureCounts.txt Human_counts/$1_KeptEBER.FeatureCounts_Human.txt > NewCounts/$1_HV.count

cat NewCounts/$1_Merge_unassigned_FeatureCounts.txt Human_counts/$1_count_FeatureCounts > NewCounts/$1_HV.count

echo "Done"
