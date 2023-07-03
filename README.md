# EBV_RNASeqPipeline

Contains the scripts for the EBV master pipeline that runs on *mandalore*.

## Quality filtering 

First we perform the quality filtering.

**How to run**

```

qsub runTrimGalore.sh Read1.fastq.gz Read2.fastq.gz

```

**Dependencies**

* TrimGalore/0.4.0
* fastqc/0.11.2
* cutadapt/1.9

## Alignment

The quality filtered reads are aligned towards the Human (Grch38) and EBV (NC_007605.1) reference genomes using STAR. We are aligning to the two reference genomes seperately to make sure we are not missing possible human-like regions within the viral reference. Afterwards dont forget to create the index for the bamfiles using samtools!

**Modification in script**:

* Paths to reference genomes
* Path to STAR executable

**How to run**

```

qsub runStar_Human_ForEBVPipeline_newSTAR.sh sample_R1_trimmed.fastq.gz sample_R2_trimmed.fastq.gz
qsub runSamtoolsindex.sh sample_NC_007605.1Aligned.sortedByCoord.out.bam

qsub runStar_Virus_ForEBVPipeline_newSTAR.sh sample_R1_trimmed.fastq.gz sample_R2_trimmed.fastq.gz
qsub runSamtoolsindex.sh sample_Grch38Aligned.sortedByCoord.out.bam

```

**Dependencies**

* STAR/2.7.2b

You will need to add Samflags to be able to filter the bamfile, doing this with samtools.

**Modification in script**:

* Path to reference genome fasta

```

qsub runAddSamflags_Grch38.sh sample_NC_007605.1Aligned.sortedByCoord.out.bam
qsub runSamtoolsindex.sh sample_Grch38Aligned.sortedByCoord.out_AddFlags.bam

qsub runAddSamflags_NC_007605.1.sh sample.bam
qsub runSamtoolsindex.sh sample_NC_007605.1Aligned.sortedByCoord.out_AddFlags.bam

```

**Dependencies**

* samtools/1.9

## EBV pipeline

The EBV pipeline filters the alignment files and outputs alignment statistics and counts. Here we filter the alignments to contain max 3 missmatches, an alignment length of more then 40 and multimapping of 10 or less. It runs in python2.

*Important*

* For first time use you need to unzip the file Homo_sapiens.GRCh38.90.sorted.ProteinCoding.gtf.gz. 

**Modification in script**:

Within runEBV_Master_pipeline_KeptEBER.sh

* The EBV Alignmnment
* The Human Alignment 
* The filtered Fastq files
* Important, you also need to change the strand information -S, is the library unstranded (0), stranded (1) or stranded reverse (2)

The bash script points to the python script EBV_Master_pipeline3.py. 

Within EBV_Master_pipeline3.py you need to change: 

* if the data is single stranded change (line 46: we now take *2 as we assume PE reads)
* line 157: featureCounts (version 2.0.0)
* line 202: featureCounts (version 2.0.0)

OBS: FeatureCount version is important here! After Version 2.0.0 the tool switched on how to handle the pair information, change the parameters for this if featurecounts is updated! 

**How to run**

```

qsub runEBV_Master_pipeline_KeptEBER.sh sample

```

**Dependencies**

Python2:
* argparse
* sys
* gzip 
* biopython
* pysam
* os

* samtools/1.9
* awk
* igvtools
* featureCounts/2.0.0
* sed
* cat


## More detailed counting after the EBV pipeline

Why? Some reads are lost across the RPMS1 due to the structure of the viral genome. Therefore we recount the EBV counts, merge to the already counted human reads again. The counting is in two steps.


### For stranded reverse sequencing

**Modification in script**:


* change path to FirstAnno (NC_007605.1-20210315_NC_changedforstrandedsequencing.tsv)
* change path to SecondAnno (NC_007605.1-20210304_RPMS1unassigned.tsv)
* change the path to your featurecounts (subread-2.0.0)
* Change path to MergeCountsFromTwoFeatureCountsRuns.py Script
* OBS! The script assumes it is a stranded reverse library, change the 2 to 1 if it is just stranded sequencing

**How to run**

```
qsub runFeatureCounts_unassignedForEBVRPMS12_StrandedRev.sh sample

```

**Dependencies**

* samtools/1.9
* featureCounts/2.0.0
* python2 

### For unstranded sequencing

**Modification in script**: 

* change path to FirstAnno (NC_007605.1-20210315_NC_changedfornon-strandedsequencing.tsv)
* change path to SecondAnno (NC_007605.1-20210304_RPMS1unassigned.tsv)
* change the path to your featurecounts (subread-2.0.0)
* Change path to MergeCountsFromTwoFeatureCountsRuns.py Script

**How to run**

```
qsub runFeatureCounts_unassignedForEBVRPMS12.sh sample

```

**Dependencies**

* samtools/1.9
* featureCounts/2.0.0
* python2



## Calc TPMS

You can calculate the tpms using the script tpm_rpkm2.R (from Andy Saurin) andrew.saurin@univ-amu.fr

**How to run**

```
Rscript tpm_rpkm2.R sample_HV.count

```

**Dependencies**
* R

## Merge FeatureCounts Table

We can use the following python script to merge the featurecounts tables from the more detailed counting

**Modification in script**:

* In case you are merging tpm tables remove the astype(int) part in pandas merge

**How to run**

```

python MergeFeatureCountTables.py *_HV.count

```

**Dependencies**
* python2
* pandas
* os
* sys