#!/usr/bin/env bash
#  __           _       _      ___ ____                      _             _       ___      _ _ _                 ___ _            _ _            
# / _\ ___ _ __(_)_ __ | |_   / _ \___ \ _  /\   /\__ _ _ __(_) __ _ _ __ | |_    / __\__ _| | (_)_ __   __ _    / _ (_)_ __   ___| (_)_ __   ___ 
# \ \ / __| '__| | '_ \| __| | | | |__) (_) \ \ / / _` | '__| |/ _` | '_ \| __|  / /  / _` | | | | '_ \ / _` |  / /_)/ | '_ \ / _ \ | | '_ \ / _ \
# _\ \ (__| |  | | |_) | |_  | |_| / __/ _   \ V / (_| | |  | | (_| | | | | |_  / /__| (_| | | | | | | | (_| | / ___/| | |_) |  __/ | | | | |  __/
# \__/\___|_|  |_| .__/ \__|  \___/_____(_)   \_/ \__,_|_|  |_|\__,_|_| |_|\__| \____/\__,_|_|_|_|_| |_|\__, | \/    |_| .__/ \___|_|_|_| |_|\___|
#                |_|                                                                                    |___/          |_|                        

# Script: 02_variant_calling_pipeline.sh
# Author: Enrique Blanco Carmona
# Year: 2019
# Project: MSc Bioinformatics UAB final thesis. 
# Objective: This pipeline is composed by many different operations: 
	# - Define variables.
	# - Generate the indexing for the reference genome.
	# - Map the reads to the reference genome to generate SAM files.
	# - Transform SAM to BAM
	# - Sort the resulting BAM files.
	# - Mark duplicates of the sorted BAM files.
	# - Join all SAMs into a master BAM.
	# - Index the master BAM file.
	# - Perform the variant calling using GATK to obtain a VCF file.
	# - Retrieve a TSV file from the VCF file containing CHROM, POS, REF, ALT, QUAL columns. 

# Define Variables.
READPATH="./data/raw_fq"
R11=$READPATH/r11.fq
R12=$READPATH/r12.fq
R21=$READPATH/r21.fq
R22=$READPATH/r22.fq
R31=$READPATH/r31.fq
R32=$READPATH/r32.fq
R41=$READPATH/r41.fq
R42=$READPATH/r42.fq
REFPATH="./data/reference"
SAMPATH="./data/sam"
BAMPATH="./data/bam"
SORTEDPATH="./data/sorted"
DUPLICATEPATH="./data/duplicate"
MASTERPATH="./data/masterbam"
VCFPATH="./data/vcf"
TSVPATH="./data/tsv"
TMP="./data/tmp"
REFGENOMENAME=chr1.fasta
GATK="./bin/gatk-4.1.3.0/gatk"
JAVA="./bin/jre1.8.0_221/bin/java"
PICARD="$JAVA -jar $./bin/picardtools/picard.jar"

# Generate the indexing for the reference genome.
bwa index $REFPATH/$REFGENOMENAME
gem-indexer $REFPATH/$REFGENOMENAME
samtools faidx $REFPATH/$REFGENOMENAME
$PICARD CreateSequenceDictionary R=$REFPATH/REFGENOMENAME O=$REFPATH/$(echo $REFGENOMENAME | cut -d . -f 1).dict

# Maps the reads onto the referenge genome to generate SAM files.
bwa mem -t 4 -o $SAMPATH/bwa/1.sam $REFPATH/$REFGENOMENAME $R11 $R12 &
bwa mem -t 4 -o $SAMPATH/bwa/2.sam $REFPATH/$REFGENOMENAME $R21 $R22 &
bwa mem -t 4 -o $SAMPATH/bwa/3.sam $REFPATH/$REFGENOMENAME $R31 $R32 &
bwa mem -t 4 -o $SAMPATH/bwa/4.sam $REFPATH/$REFGENOMENAME $R41 $R42 &
wait 

gem-mapper -t 4 -p -o $SAMPATH/gem/1.sam -I $REFPATH/$(echo $REFGENOMENAME | cut -d . -f 1).gem -1 $R11 -2 $R12 &
gem-mapper -t 4 -p -o $SAMPATH/gem/2.sam -I $REFPATH/$(echo $REFGENOMENAME | cut -d . -f 1).gem -1 $R21 -2 $R22 &
gem-mapper -t 4 -p -o $SAMPATH/gem/3.sam -I $REFPATH/$(echo $REFGENOMENAME | cut -d . -f 1).gem -1 $R31 -2 $R32 &
gem-mapper -t 4 -p -o $SAMPATH/gem/4.sam -I $REFPATH/$(echo $REFGENOMENAME | cut -d . -f 1).gem -1 $R41 -2 $R42 &
wait 

# Transform SAM to BAM
for sam in $SAMPATH/bwa/*
do
	FILENAME=(echo $sam | rev | cut -d / -f 1 | rev | cut -d . -f 1)
	$PICARD  SamFormatConverter I=$sam O=$BAMPATH/bwa/$FILENAME.bam &
done
wait

for sam in $SAMPATH/gem/*
do
	FILENAME=(echo $sam | rev | cut -d / -f 1 | rev | cut -d . -f 1)
	$PICARD  SamFormatConverter I=$sam O=$BAMPATH/gem/$FILENAME.bam &
done
wait

# Sort BAMs.

ulimit -c unlimited # Prevents errors.

for sam in $BAMPATH/bwa/*
do
	FILENAME=(echo $sam | rev | cut -d / -f 1 | rev)
	$PICARD SortSam I=$sam O=$SORTEDPATH/bwa/sorted_$FILENAME SORT_ORDER=coordinate TMP_DIR=$TMP &
done
wait 

for sam in $BAMPATH/gem/*
do
	FILENAME=(echo $sam | rev | cut -d / -f 1 | rev)
	$PICARD SortSam I=$sam O=$SORTEDPATH/gem/sorted_$FILENAME SORT_ORDER=coordinate TMP_DIR=$TMP &
done
wait 

# Mark Duplicates.

for sam in $SORTEDPATH/bwa/*
do 
	FILENAME=(echo $sam | rev | cut -d / -f 1 | rev)
	$PICARD MarkDuplicates I=$sam O=$DUPLICATEPATH/bwa/dup_$FILENAME M=$DUPLICATEPATH/bwa/metrics_$FILENAME ASO=coordinate TMP_DIR=$TMP &
done
wait

for sam in $SORTEDPATH/gem/*
do 
	FILENAME=(echo $sam | rev | cut -d / -f 1 | rev)
	$PICARD MarkDuplicates I=$sam O=$DUPLICATEPATH/gem/dup_$FILENAME M=$DUPLICATEPATH/gem/metrics_$FILENAME ASO=coordinate TMP_DIR=$TMP &
done
wait

# Generate masterbam.

samtools merge --threads 16 $MASTERPATH/bwa/masterbwa.bam $DUPLICATEPATH/bwa/sorted* &
wait
samtools merge --threads 16 $MASTERPATH/gem/mastergem.bam $DUPLICATEPATH/bwa/sorted* &
wait

# Index masterbam.

$PICARD BuildBamIndex I=$MASTERPATH/bwa/masterbwa.bam O=$MASTERPATH/bwa/masterbwa.bai &
$PICARD BuildBamIndex I=$MASTERPATH/gem/mastergem.bam O=$MASTERPATH/gem/mastergem.bai &
wait

# Perform Variant Calling
gatk HaplotypeCaller -I $MASTERPATH/bwa/masterbwa.bam -O $VCFPATH/bwa/masterbwa.vcf -R $REFPATH/$REFGENOMENAME &
gatk Haplotypecaller -I $MASTERPATH/gem/mastergem.bam -O $VCFPATH/gem/mastergem.vcf -R $REFPATH/$REFGENOMENAME &
wait 

# Obtain the TSV file to analylze.
gatk VariantsToTable -V $VCFPATH/bwa/masterbwa.vcf -O $TSVPATH/bwa/masterbwa.tsv -F CHROM -F POS -F REF -F ALT -F QUAL &
gatk VariantsToTable -V $VCFPATH/gem/mastergem.vcf -O $TSVPATH/gem/masterbwa.tsv -F CHROM -F POS -F REF -F ALT -F QUAL &
wait