#!/bin/bash

# This is a basic shell script.  It is simple a series
# of commands that are excecuted by the shell.

# Any line that begins with '#' is a comment that is ignored
# by the shell.

# The first command is to change to our working directory
# so the script can find all the files it expects
cd ~/dc_workshop/variant_calling

# Assign the name/location of our reference genome
# to a variable ($genome)
genome=data/ref_genome/ecoli_rel606.fasta

# bwa and fa index the genome (only need to do this once)
# bwa and samtools are programs that are pre-installed on our
# server
bwa index $genome
samtools faidx $genome

# Create output paths for various intermediate and result files
# The -p option means mkdir will create the whole path if it
# does not exist and refrain from complaining if it does exist
mkdir -p results/sai
mkdir -p results/sam
mkdir -p results/bam
mkdir -p results/bcf
mkdir -p results/vcf

# We will now use a loop to run the variant calling work flow
# of each of our fastq files
for fq in data/trimmed_fastq/*.fastq
do
    # We use echo for debugging/reporting to the screen
    echo "working with file $fq"

    # This command will extract the base name of the file
    # (without the path and .fastq extension) and assign it
    # to the $base variable
    base=$(basename $fq .fastq)

    echo "base name is $base"

    # We will assign various file names to variables both
    # for convenience but also to make it easier to see
    # what is going on in the sommand below
    fq=data/trimmed_fastq/$base\.fastq
    sai=results/sai/$base\_aligned.sai
    sam=results/sam/$base\_aligned.sam
    bam=results/bam/$base\_aligned.bam
    sorted_bam=results/bam/$base\_aligned_sorted.bam
    raw_bcf=results/bcf/$base\_raw.bcf
    variants=results/bcf/$base\_variants.bcf
    final_variants=results/vcf/$base\_final_variants.vcf

    # Our data are now staged.  The series of command below will run
    # the steps of the analytical workflow

    # Align the reads to the reference genome
    bwa aln $genome $fq > $sai

    # Convert the output to the SAM format
    bwa samse $genome $sai $fq > $sam

    # Convert the SAM file to BAM format
    samtools view -S -b $sam > $bam

    # Sort the BAM file
    samtools sort -f $bam $sorted_bam

    # Index the BAM file for display purposes
    samtools index $sorted_bam

    # Do the first pass on variant calling by counting
    # read coverage
    samtools mpileup -g -f $genome $sorted_bam > $raw_bcf

    # Do the SNP calling with bcftools
    bcftools view -bvcg $raw_bcf > $variants

    # Filter the SNPs for the final output
    bcftools view $variants | /usr/share/samtools/vcfutils.pl varFilter - > $final_variants

done
# and that is the end of our for loop
