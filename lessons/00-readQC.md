# Lesson QC of Sequencing read data

##Assessment and Objectives
###Declarative knowledge
*Learners will archive the following objectives in this lesson: *<br>
1. Know FASTQC is used to determine the quality of FASTQ data <br>
2. Know Trimmomatic can be used to trim and filter FASTQ data<br>
3. Recognize the following file types:<br>
    - FASTA<br>
    - FASTQ<br>
4. Know what a Phred score is and what it means<br>
5. Know older FASTQ files may have different quality scoring<br>
6. Understand a FASTQC report and what QC checks must pass in order to work with the data<br>

###Attitudinal / disposition
*Learner should conclude the workshop with the following impressions, attitudes*<br>
1. Feel confident that a data set has been adequately  QC'd<br>
2. Feel that it is important to keep notes on QC steps<br>
3.  Feels confident about how much detail to include in their documentation of work<br>

###Skills
*Learner should conclude the workshop able to do the following*<br>
1. Preview the first few lines of a FASTQ file using head or (zcat | head)<br>
2. Determine if a FASTQ file is paired or unpaired<br>
3. Determine if a FASTQ file is barcoded<br>
4. Use FASTQC to generate a report on one or multiple FASTQ files<br>
5. Interpret a FASTQC report<br>
6. Use information from the FASTQC report to decide if data is viable<br>
7. Use information from the FASTQC report to select trimming and filtering settings on Trimmommatic <br>
8. Trim and filter FASTQ files using Trimmommatic <br>
9. Verify Trimmommatic output using FASTQC<br>


## The data
We have received Fastq data files for Illumina paired end reads.  The first thing we want to do is to assess the quality of the read data and determine whether some clean-up is required.

## FastQC
FastQC (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) provides a simple way to do some quality control checks on raw sequence data coming from high throughput sequencing pipelines. It provides a modular set of analyses which you can use to give a quick impression of whether your data has any problems of which you should be aware before doing any further analysis.

The main functions of FastQC are
* Import of data from BAM, SAM or FastQ files (any variant)
* Providing a quick overview to tell you in which areas there may be problems
* Summary graphs and tables to quickly assess your data
* Export of results to an HTML based permanent report
* Offline operation to allow automated generation of reports without running the interactive application


[Sheldon to document critical outputs before tomorrow]

## Summary of FastQC Outputs
### Per Base Sequence Quality

-image placeholder example of good data-
-image placeholder example of bad data-

For each position a BoxWhisker type plot is drawn. The elements of the plot are as follows:
The central red line is the median value
The yellow box represents the inter-quartile range (25-75%)
The upper and lower whiskers represent the 10% and 90% points
The blue line represents the mean quality

The y-axis on the graph shows the quality scores. The higher the score the better the base call. The background of the graph divides the y axis into very good quality calls (green), calls of reasonable quality (orange), and calls of poor quality (red). The quality of calls on most platforms will degrade as the run progresses, so it is common to see base calls falling into the orange area towards the end of a read.
It should be mentioned that there are number of different ways to encode a quality score in a FastQ file. FastQC attempts to automatically determine which encoding method was used, but in some very limited datasets it is possible that it will guess this incorrectly (ironically only when your data is universally very good!). The title of the graph will describe the encoding FastQC thinks your file used.

## Quality Trimming and Filtering
Once we have an idea of the quality of our raw data, it is time to trim away adapters and filter out poor quality score reads.  To accomplish this task we will use Trimmomatic (http://www.usadellab.org/cms/?page=trimmomatic). 

## Vocabulary
* Phred score
* kmer
* contig
* scaffold
* N50
* adapter
* assembly
* alignment
* annotation
* paired-end
* QC
# Amazon instance ami-3c1c3454 named dataCgen-qc is set up and publicly available. It has FastQC, Trimmomatic, and sample data for the quality control tutorial. You can search for it in the public AMIs.

It is pretty trivial, but all the commands to recreate the AMI are here:
# code for AWS setup and install
mkdir data
cd data

# Fastqc
wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.2.zip
unzip fastqc_v0.11.2.zip
chmod +x FastQC/fastqc
sudo cp fastqc /usr/bin/
#Trimmomatic
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.32.zip
unzip Trimmomatic-0.32.zip
#Acquire data
mkdir SRR
cd SRR
wget apollo.huck.psu.edu/data/SRR.tar.gz
tar -zxvf SRR.tar.gz
#qc data
../FastQC/fastqc SRR447649_1.fastq # runs in < 1 min for each .fastq
java ../Trimmomatic-0.32/trimmomatic-0.32.jar SE SRR447649_1.fastq SLIDINGWINDOW:4:20 MINLEN:36 fastqc good.fq





