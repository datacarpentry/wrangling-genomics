#!/bin/bash

# run FASTQC




fastqc *.fastq
#mv zipfiles to new directory before unzipping; fastqc_reports

For loop:
for directory in */; do head -n2 $directory/summary.txt; done

AND/OR

for directory in */
> do
> cat $directory/summary.txt >>full_report.txt
> done

#grep on full report


?? cat fastq summaries into master file? 


Trimmomatic

for file in *.fastq
> do
> java -jar ~/fastq_files/Trimmomatic-0.33/trimmomatic-0.33.jar SE $file ${file}_trim.fastq SLIDINGWINDOW:4:20 MINLEN:20
> done
#do this on one and then do shell script to finish them all
First, we will create a working directory for our analysis

    mkdir dc_workshop

It will have three subdirectories

    mkdir dc_workshop/data
    mkdir dc_workshop/docs
    mkdir dc_workshop/results

The sample data we will be working with is in a hidden directory, we will
move it to our working directory

    mv .dc_sampledata_lite/untrimmed_fastq/ ~/dc_workshop/data/

Now we can get to work.

    dc_workshop/data/
    cd untrimmed_fastq/

To run the fastqc program, we call it from its location in '~/FastQC'

    ~/FastQC/fastqc *.fastq

Now, let's create a home for our results

    mkdir /home/dcuser/results/fastqc_untrimmed_fastq
    mkdir /home/dcuser/dc_workshop/fastqc_untrimmed_reads

...and move them there
    mv *.zip ~/dc_workshop/results/fastqc_untrimmed_reads/
    mv *.html ~/dc_workshop/results/fastqc_untrimmed_reads/
