# # __           _       _      ___  _  _                _      ___  ____                       _                   _    __ _ _ _   
# / _\ ___ _ __(_)_ __ | |_   / _ \| || | _    /\/\    /_\    / _ \/___ \  _ __ ___   ___   __| |   __ _ _ __   __| |  / _(_) | |_ 
# \ \ / __| '__| | '_ \| __| | | | | || |(_)  /    \  //_\\  / /_)//  / / | '_ ` _ \ / _ \ / _` |  / _` | '_ \ / _` | | |_| | | __|
# _\ \ (__| |  | | |_) | |_  | |_| |__   _|  / /\/\ \/  _  \/ ___/ \_/ /  | | | | | | (_) | (_| | | (_| | | | | (_| | |  _| | | |_ 
# \__/\___|_|  |_| .__/ \__|  \___/   |_|(_) \/    \/\_/ \_/\/   \___,_\  |_| |_| |_|\___/ \__,_|  \__,_|_| |_|\__,_| |_| |_|_|\__|
#                |_|                                                                                                               

# Script: 04_mapq_mod_and_filt.sh
# Author: Enrique Blanco Carmona
# Year: 2019
# Project: MSc Bioinformatics UAB final thesis. 
# Objective: Starting from a SAM file, all MAPQ scores are modified to three categorical thresholds: 20, 40, 60. Then, 3 files out of the original SAM 
# are generated. 


MODIFIEDPATH="./data/modified"
FILTEREDPATH="./data/filtered"
BWAPATH="./data/sam/bwa"
GEMPATH="./data/sam/gem"
GATK="./bin/gatk-4.1.3.0/gatk"
JAVA="./bin/jre1.8.0_221/bin/java"
PICARD="$JAVA -jar $./bin/picardtools/picard.jar"
TMP="./data/tmp"
REF="./data/reference/chr1.fasta"

ulimit -c unlimited

# Modify MAPQ scores BWA files. 
for file in $BWAPATH/*
do
	FILENAME=$(echo $file | rev | cut -d / -f 1 | rev)
	awk 'BEGIN {OFS = "\t"} /^@/ {print $0} !/^@/ && ($3 != "*" && $6 != "*"){$5 = 20; print $0}' > $MODIFIEDPATH/bwa/20/mod20_$FILENAME &
	awk 'BEGIN {OFS = "\t"} /^@/ {print $0} !/^@/ && ($3 != "*" && $6 != "*"){$5 = 40; print $0}' > $MODIFIEDPATH/bwa/40/mod40_$FILENAME &
	awk 'BEGIN {OFS = "\t"} /^@/ {print $0} !/^@/ && ($3 != "*" && $6 != "*"){$5 = 60; print $0}' > $MODIFIEDPATH/bwa/60/mod60_$FILENAME &
done

# Modify MAPQ scores GEM files. 
for file in $filtPATH/*
do
	FILENAME=$(echo $file | rev | cut -d / -f 1 | rev)
	awk 'BEGIN {OFS = "\t"} /^@/ {print $0} !/^@/ && ($3 != "*" && $6 != "*"){$5 = 20; print $0}' > $MODIFIEDPATH/gem/20/mod20_$FILENAME &
	awk 'BEGIN {OFS = "\t"} /^@/ {print $0} !/^@/ && ($3 != "*" && $6 != "*"){$5 = 40; print $0}' > $MODIFIEDPATH/gem/40/mod40_$FILENAME &
	awk 'BEGIN {OFS = "\t"} /^@/ {print $0} !/^@/ && ($3 != "*" && $6 != "*"){$5 = 60; print $0}' > $MODIFIEDPATH/gem/60/mod60_$FILENAME &
done

# Filter MAPQ scores BWA files.
for file in $BWAPATH/*
do
	FILENAME=$(echo $file | rev | cut -d / -f 1 | rev)
	awk 'BEGIN {OFS = "\t"} /^@/ {print $0} !/^@/ && ($3 != "*" && $6 != "*" && $5 <= 20){print $0}' > $FILTEREDPATH/bwa/20/filt20_$FILENAME &
	awk 'BEGIN {OFS = "\t"} /^@/ {print $0} !/^@/ && ($3 != "*" && $6 != "*" && $5 <= 40){print $0}' > $FILTEREDPATH/bwa/40/filt40_$FILENAME &
	awk 'BEGIN {OFS = "\t"} /^@/ {print $0} !/^@/ && ($3 != "*" && $6 != "*" && $5 <= 60){print $0}' > $FILTEREDPATH/bwa/60/filt60_$FILENAME &
done

# Filter MAPQ scores GEM files.
for file in $filtPATH/*
do
	FILENAME=$(echo $file | rev | cut -d / -f 1 | rev)
	awk 'BEGIN {OFS = "\t"} /^@/ {print $0} !/^@/ && ($3 != "*" && $6 != "*" && $5 <= 20){print $0}' > $FILTEREDPATH/gem/20/filt20_$FILENAME &
	awk 'BEGIN {OFS = "\t"} /^@/ {print $0} !/^@/ && ($3 != "*" && $6 != "*" && $5 <= 40){print $0}' > $FILTEREDPATH/gem/40/filt40_$FILENAME &
	awk 'BEGIN {OFS = "\t"} /^@/ {print $0} !/^@/ && ($3 != "*" && $6 != "*" && $5 <= 60){print $0}' > $FILTEREDPATH/gem/60/filt60_$FILENAME &
done

# TREAT MODIFIED FILES

# BWA files.
for mapq in 20 40 60
do
	for file in $MODIFIEDPATH/bwa/$mapq/mod*
	do
		FILENAME=$(echo $file | rev | cut -d / -f 1 | rev | cut -d . -f 1)
		OUTPUTPATH=$MODIFIEDPATH/$mapq
		$PICARD SamFormatConverter I=$file O=$OUTPUTPATH/$FILENAME.bam
		$PICARD SortSam I=$OUTPUTPATH/$FILENAME.bam O=$OUTPUTPATH/sorted_$FILENAME.bam SORT_ORDER=coordinate TMP_DIR=$TMP
		$PICARD MarkDuplicates I=$OUTPUTPATH/sorted_$FILENAME.bam O=$OUTPUTPATH/dup_sorted_$FILENAME.bam M=$OUTPUTPATH/metrics_sorted_$filename.bam ASO=coordinate TMP_DIR=$TMP
	done

	samtools merge --threads 16 $MODIFIEDPATH/bwa/$mapq/masterbwa_$mapq.bam $MODIFIEDPATH/bwa/$mapq/dup* 
	$PICARD BuildBamIndex I=$MODIFIEDPATH/bwa/$mapq/masterbwa_$mapq.bam O=$MODIFIEDPATH/bwa/$mapq/masterbwa_$mapq.bai
	gatk HaplotypeCaller -I $MODIFIEDPATH/bwa/$mapq/masterbwa_$mapq.bam -O $MODIFIEDPATH/bwa/vcf/masterbwa_$mapq.vcf -R $REF
	gatk VariantsToTable -V $MODIFIEDPATH/bwa/vcf/masterbwa_$mapq.vcf -O $MODIFIEDPATH/bwa/tsv/masterbwa_$mapq.tsv -F CHROM -F POS -F REF -F ALT -F QUAL
done

# GEM files.
for mapq in 20 40 60
do
	for file in $MODIFIEDPATH/gem/$mapq/mod*
	do
		FILENAME=$(echo $file | rev | cut -d / -f 1 | rev | cut -d . -f 1)
		OUTPUTPATH=$MODIFIEDPATH/$mapq
		$PICARD SamFormatConverter I=$file O=$OUTPUTPATH/$FILENAME.bam
		$PICARD SortSam I=$OUTPUTPATH/$FILENAME.bam O=$OUTPUTPATH/sorted_$FILENAME.bam SORT_ORDER=coordinate TMP_DIR=$TMP
		$PICARD MarkDuplicates I=$OUTPUTPATH/sorted_$FILENAME.bam O=$OUTPUTPATH/dup_sorted_$FILENAME.bam M=$OUTPUTPATH/metrics_sorted_$filename.bam ASO=coordinate TMP_DIR=$TMP
	done

	samtools merge --threads 16 $MODIFIEDPATH/gem/$mapq/masterbwa_$mapq.bam $MODIFIEDPATH/gem/$mapq/dup* 
	$PICARD BuildBamIndex I=$MODIFIEDPATH/gem/$mapq/masterbwa_$mapq.bam O=$MODIFIEDPATH/gem/$mapq/masterbwa_$mapq.bai
	gatk HaplotypeCaller -I $MODIFIEDPATH/gem/$mapq/masterbwa_$mapq.bam -O $MODIFIEDPATH/gem/vcf/masterbwa_$mapq.vcf -R $REF
	gatk VariantsToTable -V $MODIFIEDPATH/gem/vcf/masterbwa_$mapq.vcf -O $MODIFIEDPATH/gem/tsv/masterbwa_$mapq.tsv -F CHROM -F POS -F REF -F ALT -F QUAL
done

# TREAT FILTERED FILES.

# BWA files.
for mapq in 20 40 60
do
	for file in $FILTEREDPATH/bwa/$mapq/filt*
	do
		FILENAME=$(echo $file | rev | cut -d / -f 1 | rev | cut -d . -f 1)
		OUTPUTPATH=$FILTEREDPATH/$mapq
		$PICARD SamFormatConverter I=$file O=$OUTPUTPATH/$FILENAME.bam
		$PICARD SortSam I=$OUTPUTPATH/$FILENAME.bam O=$OUTPUTPATH/sorted_$FILENAME.bam SORT_ORDER=coordinate TMP_DIR=$TMP
		$PICARD MarkDuplicates I=$OUTPUTPATH/sorted_$FILENAME.bam O=$OUTPUTPATH/dup_sorted_$FILENAME.bam M=$OUTPUTPATH/metrics_sorted_$filename.bam ASO=coordinate TMP_DIR=$TMP
	done

	samtools merge --threads 16 $FILTEREDPATH/bwa/$mapq/masterbwa_$mapq.bam $FILTEREDPATH/bwa/$mapq/dup* 
	$PICARD BuildBamIndex I=$FILTEREDPATH/bwa/$mapq/masterbwa_$mapq.bam O=$FILTEREDPATH/bwa/$mapq/masterbwa_$mapq.bai
	gatk HaplotypeCaller -I $FILTEREDPATH/bwa/$mapq/masterbwa_$mapq.bam -O $FILTEREDPATH/bwa/vcf/masterbwa_$mapq.vcf -R $REF
	gatk VariantsToTable -V $FILTEREDPATH/bwa/vcf/masterbwa_$mapq.vcf -O $FILTEREDPATH/bwa/tsv/masterbwa_$mapq.tsv -F CHROM -F POS -F REF -F ALT -F QUAL
done

# GEM files.
for mapq in 20 40 60
do
	for file in $FILTEREDPATH/gem/$mapq/filt*
	do
		FILENAME=$(echo $file | rev | cut -d / -f 1 | rev | cut -d . -f 1)
		OUTPUTPATH=$FILTEREDPATH/$mapq
		$PICARD SamFormatConverter I=$file O=$OUTPUTPATH/$FILENAME.bam
		$PICARD SortSam I=$OUTPUTPATH/$FILENAME.bam O=$OUTPUTPATH/sorted_$FILENAME.bam SORT_ORDER=coordinate TMP_DIR=$TMP
		$PICARD MarkDuplicates I=$OUTPUTPATH/sorted_$FILENAME.bam O=$OUTPUTPATH/dup_sorted_$FILENAME.bam M=$OUTPUTPATH/metrics_sorted_$filename.bam ASO=coordinate TMP_DIR=$TMP
	done

	samtools merge --threads 16 $FILTEREDPATH/gem/$mapq/masterbwa_$mapq.bam $FILTEREDPATH/gem/$mapq/dup* 
	$PICARD BuildBamIndex I=$FILTEREDPATH/gem/$mapq/masterbwa_$mapq.bam O=$FILTEREDPATH/gem/$mapq/masterbwa_$mapq.bai
	gatk HaplotypeCaller -I $FILTEREDPATH/gem/$mapq/masterbwa_$mapq.bam -O $FILTEREDPATH/gem/vcf/masterbwa_$mapq.vcf -R $REF
	gatk VariantsToTable -V $FILTEREDPATH/gem/vcf/masterbwa_$mapq.vcf -O $FILTEREDPATH/gem/tsv/masterbwa_$mapq.tsv -F CHROM -F POS -F REF -F ALT -F QUAL
done



