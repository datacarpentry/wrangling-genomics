---
layout: page
title: Setup
permalink: /setup/
---

## Table of Content

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

Most of the genomics lessons from data carpentry currently use amazon cloud.

We don't currently know if we'll keep using amazon cloud or not.

## Required software

### FastQC

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) provides a simple way to do some quality control checks on raw sequence data coming from high throughput sequencing pipelines. It provides a modular set of analyses which you can use to give a quick impression of whether your data has any problems of which you should be aware before doing any further analysis.

FastQC is available for Linux, MacOS and Windows.

### Trimmomatic

[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) is a java based program that can remove sequencer specific reads and nucleotides that fall below a certain threshold. Trimmomatic can be multithreaded to run quickly.

Trimmomatic is available for Linux, MacOS and Windows.

### BWA

[Bwa](https://github.com/lh3/bwa) is a software package for mapping DNA sequences against a large reference genome, such as the human genome.

Bwa is available for Linux and MacOS.

#### Installation instructions:
Reference: The Biostars Handbook

```brew install bwa```
or

```conda install -y bwa```
All other platforms:
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

[SAMtools](https://github.com/samtools/samtools) is a suite of programs for interacting with high-throughput sequencing data. Samtools can read/write/edit/index/view SAM/BAM/CRAM format.

SAMtools is available for Linux and MacOS.


### bcftools

[BCFtools](https://github.com/samtools/bcftools) is a program for variant calling and manipulating files in the Variant Call Format (VCF) and its binary counterpart BCF.

BCFtools is available for Linuc and MacOS

### IGV

[IGV](http://software.broadinstitute.org/software/igv/) is a high-performance visualization tool for interactive exploration of large, integrated genomic datasets.

IGV is available for Linux, MacOS and Windows.

## Required Data

You will also need to download a data tarball of a reference genome and fastq files for *E. coli*:

* Download [variant_calling.tar.gz](./variant_calling.tar.gz), you can also use the unix command 'wget' to obtain the copied link.
* Once downloaded unpack it with `tar xzf variant_calling.tar.gz`
