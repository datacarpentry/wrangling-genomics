# Lesson QC of Sequence Read Data

##Assessment and Objectives
###Declarative knowledge
*Learners will archive the following objectives in this lesson:*<br>
1. Know FASTQC is used to determine the quality of FASTQ data <br>
2. Know Trimmomatic can be used to trim and filter FASTQ data<br>
3. Recognize the following file types:<br>
* FASTA
* FASTQ
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

## Amazon instance 
An Amazon Linux AMI ami-3c1c3454 named dataCgen-qc is set up and publicly available. It has FastQC, Trimmomatic, and sample data for the quality control tutorial. You can search for it in the public AMIs.  It is assumed that you know how to launch an AWS instance and have the private key downloaded to your computer.

### Setup
You will need:
* A t2.medium instance of mi-3c1c3454
* The pem file from you public/private key pair
* The public IP address of your AWS instance

### Connecting
* `<pemfile.pem>` is the name of your locally saved private key file
* `<public IP>` is the IP address of your AWS instance
```
ssh -i <pemfile.pem> ec2-user@<public IP>
Last login: Wed Mar 25 10:14:10 2015 from ool-addc1c9a.static.optonline.net

       __|  __|_  )
       _|  (     /   Amazon Linux AMI
      ___|\___|___|

https://aws.amazon.com/amazon-linux-ami/2014.09-release-notes/
18 package(s) needed for security, out of 147 available
Run "sudo yum update" to apply all updates.
Amazon Linux version 2015.03 is available.
```
## The data
We have received Fastq data files for Illumina paired end reads.  The first thing we want to do is to assess the quality of the read data and determine whether some clean-up is required.  FastQC has been run on one of the FastQ files.

* The data are on the folder `~/data/SRR/`

```
ls -al ~/data/SRR
total 1523764
drwxrwxr-x 2 ec2-user ec2-user      4096 Mar 25 01:37 .
drwxrwxr-x 5 ec2-user ec2-user      4096 Mar 25 01:36 ..
-rw-r--r-- 1 ec2-user ec2-user 453225654 Jul 22  2014 SRR447649_1.fastq
-rw-rw-r-- 1 ec2-user ec2-user    370420 Mar 25 10:17 SRR447649_1_fastqc.html
-rw-rw-r-- 1 ec2-user ec2-user    675772 Mar 25 10:17 SRR447649_1_fastqc.zip
-rw-r--r-- 1 ec2-user ec2-user 192198334 Jul 22  2014 SRR447649_2.fastq
-rw-r--r-- 1 ec2-user ec2-user 453225654 Jul 22  2014 SRR447649_3.fastq
-rw-r--r-- 1 ec2-user ec2-user 230298876 Jul 22  2014 SRR519926_1.fastq
-rw-r--r-- 1 ec2-user ec2-user 230298876 Jul 22  2014 SRR519926_2.fastq
```

* The FastQC software is installed in the folder `~/data/FastQC`
```
 ls -al ~/data/FastQC/
total 820
drwxrwxr-x 8 ec2-user ec2-user   4096 Mar 25 01:30 .
drwxrwxr-x 5 ec2-user ec2-user   4096 Mar 25 01:36 ..
drwxrwxr-x 2 ec2-user ec2-user   4096 Jun  6  2014 Configuration
-rwxrwxr-x 1 ec2-user ec2-user  12802 Jun  4  2014 fastqc
-rw-rw-r-- 1 ec2-user ec2-user   2238 Mar 21  2012 fastqc_icon.ico
drwxrwxr-x 5 ec2-user ec2-user   4096 Jun  6  2014 Help
-rw-rw-r-- 1 ec2-user ec2-user   5611 May  7  2014 INSTALL.txt
-rw-rw-r-- 1 ec2-user ec2-user  50147 Feb 24  2014 jbzip2-0.9.jar
-rw-rw-r-- 1 ec2-user ec2-user  35821 Mar 21  2012 LICENSE.txt
drwxrwxr-x 3 ec2-user ec2-user   4096 Jun  6  2014 net
drwxrwxr-x 3 ec2-user ec2-user   4096 Jun  6  2014 org
-rw-rw-r-- 1 ec2-user ec2-user   2143 Mar 21  2012 README.txt
-rw-rw-r-- 1 ec2-user ec2-user  31771 Jun  6  2014 RELEASE_NOTES.txt
-rw-rw-r-- 1 ec2-user ec2-user    101 May  6  2014 run_fastqc.bat
-rw-rw-r-- 1 ec2-user ec2-user 644848 Feb 14  2014 sam-1.103.jar
drwxrwxr-x 3 ec2-user ec2-user   4096 Jun  6  2014 Templates
drwxrwxr-x 3 ec2-user ec2-user   4096 Jun  6  2014 uk
```
The file `fastqc` is the one that will be executed.  Make sure it is executable:
```
chmod +x ~/data/SRR/FastQC/fastqc
```

<!--
### Set up apache web server (for
fastQC results)  
**NEVERMIND, just have them download the files to their desktop**
* We will set up our AWS instance as a web server so we can view and share our FastQC results
* This distribution of Linux used the [redhat yum package manager](https://access.redhat.com/solutions/9934)
* It is very easy to install the [apache web server software](http://httpd.apache.org/)
* **Note: ** you need to use root-level access to install apache (SHELDON says: we should pre-install this on the AMI)

`sudo yum install httpd`

* Once the software is installed, the (document root) is `/var/www/html`, which is where you will copy content you wish to publish on the web
** **Note: ** You need root (sudo) access to copy files here (SHELDON says: we should probably just make ec2-user the owner or usr a symlink to ~ec2-user/data/SRR/fastQC or something like that)
-->

## FastQC
FastQC (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) provides a simple way to do some quality control checks on raw sequence data coming from high throughput sequencing pipelines. It provides a modular set of analyses which you can use to give a quick impression of whether your data has any problems of which you should be aware before doing any further analysis.

The main functions of FastQC are
* Import of data from BAM, SAM or FastQ files (any variant)
* Providing a quick overview to tell you in which areas there may be problems
* Summary graphs and tables to quickly assess your data
* Export of results to an HTML based permanent report
* Offline operation to allow automated generation of reports without running the interactive application

## Running FastQC
Note: this has been run already for the first file
```
cd ~/data/SRR
../data/FastQC/fastqc SRR447649_1.fastq # runs in < 1 min for each .fastq
```
## Viewing the FastQC output


<!-- Don't need this!
## Running FastQC
```
cd ~/data/SRR
../data/FastQC/fastqc SRR447649_1.fastq # runs in < 1 min for each .fastq
```

## Use web server to view the results
### Put the data where the web server can see it
placeholder more helpful words here
```
sudo cp SRR447649_1_fastqc* /var/www/html
cd /var/www/html
sudo unzip SRR447649_1_fastqc.zip
```

### From the web browser
On your (local) computer, point any web browser to the address below,
where `<public IP>` is the ip address of your AWS instance.

`http://<public IP>/SRR447649_1_fastqc.html`

image placeholder fastqc front page SUCCESS!
-->

## Summary of FastQC Outputs
### Per Base Sequence Quality
image placeholder example of good data-<br>
image placeholder example of bad data-
* what QC steps can be taken based on this report
### more
* what QC steps can be taken based on this report
### still more
* what QC steps can be taken based on this report

# Placeholder splitting paired read files
* Do we need to?
* How do we do it?


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



# (Optional) steps for software installation on the AWS instance
```
mkdir data
cd data
```

## FastQC
```
wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.2.zip
unzip fastqc_v0.11.2.zip
chmod +x FastQC/fastqc
sudo cp fastqc /usr/bin/
```
## Trimmomatic
    wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.32.zip
    unzip Trimmomatic-0.32.zip
    
## Acquiring data
    mkdir SRR
    cd SRR
    wget apollo.huck.psu.edu/data/SRR.tar.gz
    tar -zxvf SRR.tar.gz
    
## Trimming reads
    java ../Trimmomatic-0.32/trimmomatic-0.32.jar SE SRR447649_1.fastq SLIDINGWINDOW:4:20 MINLEN:36 fastqc good.fq







