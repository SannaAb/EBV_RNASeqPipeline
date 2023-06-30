#!/usr/bin/py

# This is the entire pipeline that takes the input of a Alignment File, Filters the alignmenfile in the userdefined parameters. Runs the statistics by inputting the belonging FastqFile and then also run featureCount for the Entire alignment


# Obs, the scripts works in the environment Sanna2 which contains pysam and biopython

import argparse
import sys
import gzip
from Bio import SeqIO
import pysam
import os


def parseArgs():
    parser = argparse.ArgumentParser(description='Collects the EBV ppm, counts for human and EBV together with a filtering step')
    parser.add_argument('-IH', dest='Hbaminput',help='Human Bamfile containing alignment performed by a splice aware aligner, need index in the same folder (required)',required=True)
    parser.add_argument('-IE', dest='Ebaminput',help='EBV Bamfile containing alignment performed by a splice aware aligner, need index in the same folder (required)',required=True)
    parser.add_argument('-O', dest='OutputFile', help='Outputfilename the prefix of your outputs together with the output folder (required)', required=True)
    parser.add_argument('-F', dest='FastqFile', help='Fastq File, either read 1 or read 2. OBS need to be paired end! This is used for calculating the ppm (required)', required=True)
    parser.add_argument('-MM', dest='missmatchTreshold', help='The treshold that the user want to filter the alignment based on missmatch Treshold', default=3, type=int)
    parser.add_argument('-AL', dest='alignmentLenght', help='The minimal lenght that the read aligning towards the genome have', default=40, type=int)
    parser.add_argument('-MU', dest='Multimapped', help='The maximum amount of times a read is allowed to mapp', default=10, type=int)
    parser.add_argument('--filterEBER', help='The minimal lenght that the read aligning towards the genome have', action='store_true') # If you add filter EBER then we filter eber from the alignment, otherwise skip
    parser.add_argument('-S', dest='strandness', help='Is the data unstranded (0), stranded (1) or stranded reverse (2)',required=True,  type=int)
    
    arguments=parser.parse_args(sys.argv[1:])
    return arguments

def CountTotalAmountofStartReads(FastqFile):
    """
    This part calculates the total amount of reads in the fastq file, takes both zipped and unzipped file
    """
    TotalAmountOfReads = 0
    print "Calculating amount of raw reads from the fastq..."
    if FastqFile.endswith(".gz"): # The file is zipped
        with gzip.open(FastqFile, "r") as handle:
            for record in SeqIO.parse(handle, "fastq"):
                TotalAmountOfReads += 1
    else: # The file is not zipped
        handle = open(fastqfile,"r")
        for record in SeqIO.parse(handle, "fastq"):
            TotalAmountOfReads += 1
        handle.close()
    TotalAmountOfReads = TotalAmountOfReads * 2  # The data is paired end... I assume
    return TotalAmountOfReads

def FilterEBV(Ebaminput,OutputFile,missmatchTreshold,alignmentLenght,Multimapped,TotalAmountOfReads,filterEBER):
    """
    This Part filters the Alignmentfiles from the EBV Alignment
    """
    FilteredFilename = "%s.%sMissmatch_%snt_%sMultimapped.EBV.bam" %(OutputFile,missmatchTreshold,alignmentLenght,Multimapped)

    EBERCoordsBED = os.path.abspath(__file__).split("Code/EBV_Master_pipeline3.py")[0] + "DB/NC_WithAJ50anno_EBER.coords.bed"
    #EBERCoordsBED = "/medstore/projects/B18-045/db/EBV/Annotation/NC_WithAJ50anno_EBER.coords.bed"
    pathtoNCgenome = os.path.abspath(__file__).split("Code/EBV_Master_pipeline3.py")[0] + "DB/NC_007605.1.fasta"
    #pathtoNCgenome = "/medstore/projects/B18-045/db/EBV/Sequences/NC_007605.1.fasta"
    command = "samtools view %s -h | awk \"\$1 ~ /@/" %(Ebaminput)
    missmatchpart = ""
    # Add the missmatches
    print "Filtering the alignment file: missmatch: %s, AlignmentLength: %s, Multimapping: %s" %(missmatchTreshold,alignmentLenght,Multimapped)
    for i in range(0,missmatchTreshold+1): # Create the missmatchloop, a string for filtering the missmatches
        missmatchpart = missmatchpart+"|| \$16 == \\\"NM:i:%s\\\"" %(i)
    command = command + missmatchpart + "{print \$0}\""
    # Add the length
    command = command + "| awk \"\$1 ~ /@/ ||  length(\$10) >= %s\"| awk \"\$1 ~ /@/ " %(alignmentLenght)
    # Add the multimapping
    multimappedpart=""
    for i in range(0,Multimapped+1):
        multimappedpart = multimappedpart+"|| \$0 ~\\\"NH:i:%s\\\" " %(i)
    command = command + multimappedpart + "{print \$0}\" |  samtools view -Sb > %s" % FilteredFilename
    # We probably need a check to make sure that the alignment file is not empty before or is empty after filtering
    os.system(command)
    command = "samtools index %s" % FilteredFilename
    os.system(command)
    # create the tdf
    tdfout = FilteredFilename.split(".bam")[0] + ".tdf"
    command = "igvtools count %s %s %s" %(FilteredFilename,tdfout,pathtoNCgenome)
    os.system(command)
    if filterEBER:
        print "Filtering EBER..."
        FilteredFilenameWithoutEBER = FilteredFilename.split(".bam")[0] + ".WithoutEBER.bam"
        # If we have the flag for filtering EBER we are removing EBER
        command = "samtools view -b %s -L %s -U %s" %(FilteredFilename,EBERCoordsBED,FilteredFilenameWithoutEBER)
        os.system(command)
        command = "samtools index %s" %FilteredFilenameWithoutEBER
        os.system(command)
        tdfout =  FilteredFilename.split(".bam")[0] + ".WithoutEBER.tdf"
        FilteredFilename = FilteredFilenameWithoutEBER
        command = "igvtools count %s %s %s" %(FilteredFilename,tdfout,pathtoNCgenome)
        os.system(command)
    return FilteredFilename

def CalCulateMappingStatisticsForEB(FilteredFilename,TotalAmountOfReads,OutputFile):
    """
    This part calculates the amount of mapped reads towards EBV, The ppms and the amount of mapped reads in the RPMS1 region
    """
    RPSM1Genome = "NC_007605.1"
    RPMS1Coordsstart = 138352
    RPMS1Coordsend = 160531
    samf=pysam.AlignmentFile(FilteredFilename, "rb")
    resdict = {}
    resdict["TotalAmountOfReads"]=TotalAmountOfReads
    outmappingStatistics = OutputFile + ".EBV_RPMS1.Stats.txt"
    # First we check the total amount of mapped reads in EBV
    print "Calculating amount of reads mapped to EBV"
    totalamountofmappedreads = 0
    ReadsacrossEntireGenome=samf.fetch()
    for read in ReadsacrossEntireGenome:
        if not read.is_unmapped:
            totalamountofmappedreads +=1
    resdict["TotalAmountOfMappedReadsEBV"] = totalamountofmappedreads
    # Then we check how many reads that mapps in the RPMS1 coords
    print "Calculating amount of reads mapping to RPMS1"
    totalamountofmappedreadsRPMS1 = 0
    ReadsacrossRPMS1=samf.fetch(RPSM1Genome, RPMS1Coordsstart, RPMS1Coordsend)
    for read in ReadsacrossRPMS1:
        totalamountofmappedreadsRPMS1 += 1
    resdict["TotalAmountOfMappedReadsRPMS1"] = totalamountofmappedreadsRPMS1
    samf.close()
    # Lets calculate the fraction of the reads mapping in the RPM1 region compared to the entire region
    print "Calculating fraction of RPMS1 reads of the total EBV reads"
    try: # We only have 0
        resdict["FractionInRPMS1ofEBVreads"] =resdict["TotalAmountOfMappedReadsRPMS1"] / float(resdict["TotalAmountOfMappedReadsEBV"])
    except ZeroDivisionError:
        resdict["FractionInRPMS1ofEBVreads"] = 0
    # Finally lets calculate the EBV ppm
    print "Calculating EBV PPM"
    try:
        resdict["EBVPPM"] = float(resdict["TotalAmountOfMappedReadsEBV"])/TotalAmountOfReads*1000000
    except ZeroDivisionError:
        resdict["EBVPPM"] = 0
    neworder = [2,4,3,1,0] # This part only reorders the list so the output will be nicer
    valuelist = resdict.values()
    valuelist = [valuelist[i] for i in neworder]
    keylist = resdict.keys()
    keylist = [keylist[i] for i in neworder]
    # And output It
    with open(outmappingStatistics, "w") as out:
        print >> out, "\t".join(keylist)
        print >> out, "\t".join(str(x) for x in valuelist) # This part is needed for converting the ints to strings

def CalculateCountEBV(FilteredFilename, strandness, OutputFile):
    """
    This part runs feature counts for the EBV to get the gene count for EBV
    """
    print "Counting genes in EBV"
    
    pathtoanno = os.path.abspath(__file__).split("Code/EBV_Master_pipeline3.py")[0] + "DB/NC_007605.1-201909171.gff3"

    outputFeatureCountsEBV = OutputFile + ".FeatureCounts_EBV.txt"
    if strandness not in [0,1,2]:
        print("Library needs to be either of 0,1,2. Killing it")
        sys.exit()
    else:
        command = "/home/xabras/Programs/subread-2.0.0-source/bin/featureCounts  -a %s -t exon -g gene -M -p -s %s -J -T 1 -o %s %s" %(pathtoanno,strandness,outputFeatureCountsEBV,FilteredFilename)
        os.system(command)
    return outputFeatureCountsEBV

def FilteringHuman(Hbaminput,OutputFile,missmatchTreshold,alignmentLenght,Multimapped):
    """
    This part filters the human Alignment based on the same paramters as the Viral alignment
    """
    print "Filtering the Human alignment"
    FilteredFilenameHuman = "%s.%sMissmatch_%snt_%sMultimapped.Human.bam" %(OutputFile,missmatchTreshold,alignmentLenght,Multimapped)
    command = "samtools view %s -h | awk \"\$1 ~ /@/" %(Hbaminput)
    missmatchpart = ""
    # Add the missmatches
    print "Filtering the alignment file: missmatch: %s, AlignmentLength: %s, Multimapping: %s" %(missmatchTreshold,alignmentLenght,Multimapped)
    for i in range(0,missmatchTreshold+1): # Create the missmatchloop, a string for filtering the missmatches
        missmatchpart = missmatchpart+"|| \$16 == \\\"NM:i:%s\\\"" %(i)
    command = command + missmatchpart + "{print \$0}\""
    # Add the length
    command = command + "| awk \"\$1 ~ /@/ ||  length(\$10) >= %s\"| awk \"\$1 ~ /@/ " %(alignmentLenght)
    # Add the multimapping
    multimappedpart=""
    for i in range(0,Multimapped+1):
        multimappedpart = multimappedpart+"|| \$0 ~\\\"NH:i:%s\\\" " %(i)
    command = command + multimappedpart + "{print \$0}\" |  samtools view -Sb > %s" % FilteredFilenameHuman
    # We probably need a check to make sure that the alignment file is not empty before or is empty after filtering
    os.system(command)
    command = "samtools index %s" % FilteredFilenameHuman
    os.system(command)
    return FilteredFilenameHuman

def CalculateCountHuman(FilteredFilenameHuman, strandness, OutputFile):
    """
    This part runs feature counts for the EBV to get the gene count for EBV
    """
    print "Counting genes in Human, obs we are only counting protein coding genes"

    pathtoanno = os.path.abspath(__file__).split("Code/EBV_Master_pipeline3.py")[0] + "DB/Homo_sapiens.GRCh38.90.sorted.ProteinCoding.gtf"

    #pathtoanno = "/medstore/projects/B18-045/db/HG38/Homo_sapiens.GRCh38.90.sorted.ProteinCoding.gtf"
    outputFeatureCountsHuman = OutputFile + ".FeatureCounts_Human.txt"

    if strandness not in [0,1,2]:
        print("Error: Library needs to be either of 0,1,2. Exiting now")
        sys.exit()
    else:
        command = "/home/xabras/Programs/subread-2.0.0-source/bin/featureCounts  -a %s -t exon -g gene_name -M -p -s %s -J -T 1 -o %s %s" %(pathtoanno,strandness,outputFeatureCountsHuman,FilteredFilenameHuman)
        os.system(command)
    return outputFeatureCountsHuman

def CalculateTPM(outputFeatureCountsEBV,outputFeatureCountsHuman,OutputFile):
    """
    This part calculates the TPMS by combining the counts of EBV and the protein coding human genes
    """
    print "Combining Viral Counts with human"
    NormalizedCountsScript = os.path.abspath(__file__).split("Code/EBV_Master_pipeline3.py")[0] + "Code/tpm_rpkm2.R"

    #NormalizedCountsScript = "/home/xabras/Scripts/tpm_rpkm2.R" # I just needed a slight modification for the tpm_rpkm script to keep the structure if we only substract one column which is true in this script as we run the samples one by one
    outputFeatureCountsCombined = OutputFile + ".FeatureCounts_VH.txt"
    command = "sed '1,2d' -i %s" %outputFeatureCountsHuman # This part removed the header rows from the human count matrix
    os.system(command)
    command = "cat %s %s > %s" %(outputFeatureCountsEBV,outputFeatureCountsHuman,outputFeatureCountsCombined)
    os.system(command)
    print "Calculating the TPM"
    # Now lets calculate the tpm using the Rscript ..
    command = "Rscript %s %s" %(NormalizedCountsScript,outputFeatureCountsCombined)
    os.system(command)
    tpmout = outputFeatureCountsCombined.split(".txt")[0] + "_tpm.txt"
    rpkmout = outputFeatureCountsCombined.split(".txt")[0] + "_rpkm.txt"
    return (tpmout,rpkmout)

def main(Hbaminput,Ebaminput,OutputFile,FastqFile,missmatchTreshold,alignmentLenght,Multimapped,filterEBER, strandness):
    # Make the outputFolders
    TotalAmountOfReads=CountTotalAmountofStartReads(FastqFile)
    FilteredFilename=FilterEBV(Ebaminput,OutputFile,missmatchTreshold,alignmentLenght,Multimapped,TotalAmountOfReads,filterEBER)
    CalCulateMappingStatisticsForEB(FilteredFilename,TotalAmountOfReads,OutputFile)
    outputFeatureCountsEBV=CalculateCountEBV(FilteredFilename,strandness, OutputFile)
    FilteredFilenameHuman=FilteringHuman(Hbaminput,OutputFile,missmatchTreshold,alignmentLenght,Multimapped)
    outputFeatureCountsHuman=CalculateCountHuman(FilteredFilenameHuman,strandness, OutputFile)
    (tpmout,rpkmout)=CalculateTPM(outputFeatureCountsEBV,outputFeatureCountsHuman,OutputFile)


if __name__=='__main__':
  arguments=parseArgs()
  main(arguments.Hbaminput, arguments.Ebaminput, arguments.OutputFile, arguments.FastqFile, arguments.missmatchTreshold, arguments.alignmentLenght,arguments.Multimapped, arguments.filterEBER, arguments.strandness)
