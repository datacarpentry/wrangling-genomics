---
layout: page
title: Setup
permalink: /setup/
---

## Table of Contents

* [Amazon Cloud](#amazon-cloud)
* [Required software](#required-software)
    * [FastQC](#fastqc)
    * [Trimmomatic](#trimmomatic)
    * [BWA](#bwa)
    * [SAMtools](#samtools)
    * [BCFtools](#bcftools)
    * [IGV](#igv)
* [Required data](#required-data)

## Amazon Cloud

Most of the genomics lessons from Data Carpentry currently use Amazon Cloud.

We don't currently know if we'll keep using Amazon Cloud or not.

## Required software

### FastQC

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) provides a simple way to do some quality control checks on raw sequence data coming from high throughput sequencing pipelines. It provides a modular set of analyses which you can use to give a quick impression of whether your data has any problems of which you should be aware before doing any further analysis.

FastQC is available for Linux, MacOS and Windows.

#### Install Instructions:
Reference: The Biostars Handbook

```brew install fastqc```
or

```conda install -y fastqc```

#### Source code installation

**A note on `curl` and `wget`**  
In these instructions, the command `curl` is used to download files from the internet. This command supports downloading files from FTP and HTTP(s). The `-O` paremeter ensures that the downloaded file gets saved on disk with the same filename as the original. There are other commandline tools that can also be used to download data, for example `wget` (which you should use without the `-O` flag). `wget` also supports recursive download (with the parameter `-r`), allowing you to download content from a directory or folder. Not all computers have both tools installed, though.

This is helpful when one wants to understand what type of files come with fastqc
````
cd ~/src
curl -O http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip
unzip fastqc_v0.11.5.zip

# Link the fastqc executable to the ~/bin folder that
# you have already added to the path.
ln -sf ~/src/FastQC/fastqc ~/bin/fastqc

# Due to what seems a packaging error
# the executable flag on the fastqc program is not set.
# We need to set it ourselves.
chmod +x ~/bin/fastqc
```
Test installation by running:

```fastqc -h```

### Trimmomatic

[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) is a java based program that can remove sequencer specific reads and nucleotides that fall below a certain threshold. Trimmomatic can be multithreaded to run quickly.

Trimmomatic is available for Linux, MacOS and Windows.

#### Installation Instructions:
Reference: The Biostars Handbook

```brew install trimmomatic```
or
```conda install -y trimmomatic```
#### Source Code Installation

cd ~/src
curl -O http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip
unzip Trimmomatic-0.36.zip

# The program can be invoked via
java -jar ~/src/Trimmomatic-0.36/trimmomatic-0.36.jar

# The ~/src/Trimmomatic-0.36/adapters/ directory contains
# Illumina specific adapter sequences.
ls ~/src/Trimmomatic-0.36/adapters/
```
#### How to run trimmomatic

Unfortunately running trimmomatic is as user unfriendly as it gets. To run it we "simply" type:

```java -jar ~/src/Trimmomatic-0.36/trimmomatic-0.36.jar```
That gets old very quickly. To simplify the invocation create a script in the ~/bin folder:

```
echo '#!/bin/bash' > ~/bin/trimmomatic
echo 'java -jar ~/src/Trimmomatic-0.36/trimmomatic-0.36.jar $@' >> ~/bin/trimmomatic
chmod +x ~/bin/trimmomatic
```
#### Test installation by running:

```trimmomatic```

### BWA

[Bwa](https://github.com/lh3/bwa) is a software package for mapping DNA sequences against a large reference genome, such as the human genome.

Bwa is available for Linux and MacOS.

#### Installation instructions:
Reference: The Biostars Handbook

```brew install bwa```
or

```conda install -y bwa```

##### Installation from source:
```
cd ~/src
curl -OL http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.15.tar.bz2
tar jxvf bwa-0.7.15.tar.bz2
cd bwa-0.7.15
make
export PATH=~/src/bwa-0.7.15:$PATH
```

#### Test installation:

Run:
```
bwa
```
BWA has a nicely formatted manual:

```man ~/src/bwa-0.7.12/bwa.1 ```

### SAMtools

[SAMtools](https://github.com/samtools/samtools) is a suite of programs for interacting with high-throughput sequencing data. SAMtools can read/write/edit/index/view SAM/BAM/CRAM format.

SAMtools is available for Linux and MacOS.

#### Installation Instructions
Reference: The Biostars Handbook

```brew install samtools```
or
```conda install -y samtools```
Note:
SAMtools has changed the command line invocation (for the better). But this means that most of the tutorials on the web indicate an older and obsolete usage.

Use only SAMtools 1.3 or later.

#### Source code installation

```
cd ~/src
curl -OkL https://github.com/samtools/samtools/releases/download/1.3/samtools-1.3.tar.bz2
tar jxvf samtools-1.3.tar.bz2
cd samtools-1.3
make

# Add directory to the path if necessary
echo export `PATH=~/src/samtools-1.3:$PATH` >> ~/.bashrc
source ~/.bashrc
```

#### Test that the installation succeeded:

`samtools`

SAMtools has a nicely formatted manual:
```
man ~/src/samtools-1.3/samtools.1
```

### BCFtools

[BCFtools](https://github.com/samtools/bcftools) is a program for variant calling and manipulating files in the Variant Call Format (VCF) and its binary counterpart BCF.

BCFtools is available for Linux and MacOS

#### Installation Instructions

```brew install bcftools```
or
```
conda install bcftools```

#### Install from source

```
cd ~/src
curl -OkL https://github.com/samtools/bcftools/releases/download/1.5/bcftools-1.5.tar.bz2
tar jxvf bcftools-1.5.tar.bz2
cd bcftools-1.5
make

# Add directory to the path if necessary
echo export `PATH=~/src/bcftools-1.5:$PATH` >> ~/.bashrc
source ~/.bashrc
```

#### Test that the installation succeeded:

`bcftools`

### IGV

[IGV](http://software.broadinstitute.org/software/igv/) is a high-performance visualization tool for interactive exploration of large, integrated genomic datasets.

IGV is available for Linux, MacOS and Windows.

## Required Data

You will also need to download a data tarball of a reference genome and fastq files for *E. coli*, available here  [variant_calling.tar.gz](./variant_calling.tar.gz):

```
curl -O http://www.datacarpentry.org/wrangling-genomics/variant_calling.tar.gz

```

After getting the file, a good practice is to make sure that we actually got the correct file, and it's not for example corrupted or a different version than the one originally intended. The use of **checksums** is one of the most commonly used methods for ensuring that. A checksum is a 32-character string that is the unique signature of a file. In order to facilitate this check, the providers of datasets also give a file with the checksum number so that the check can be performed. In our case, the checksum is available in the file [variant_calling.md5](./variant_calling.md5), and it contains the following line:

```
55936b1e819246236535798d8a36c134  variant_calling.tar.gz
```

The 32-character string on the left is the checksum of the file listed on the right. Let's download this file as well with wget (`curl -O http://www.datacarpentry.org/wrangling-genomics/variant_calling.md5`) and check if everything is ok with the following command (**Important note**: _both data file and checksum file should be in the same directory_):

```
md5sum -c variant_calling.md5
```

This command should print out a single line, stating `variant_calling.tar.gz: OK`. If that is not true, this means that we have downloaded either the incorrect data file or the md5 file does not correspond to the file we currently have.

Finally, once we are sure that we downloaded the correct file, we can unpack it with `tar xzf variant_calling.tar.gz`.
