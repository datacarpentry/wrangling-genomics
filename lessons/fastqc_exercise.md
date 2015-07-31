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
