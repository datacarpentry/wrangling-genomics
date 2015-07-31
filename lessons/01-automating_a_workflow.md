# Lesson Automating a variant calling workflow

Quality Control of NGS Data
===================

Learning Objectives:
-------------------
#### What's the goal for this lesson?

* Use a series of command line tools to perform a variant calling workflow
* Use a For loop from the previous lesson to help automate repetitive tasks
* Group a series of sequential commands into a script to automate a workflow

To get started with this lesson, we will need to grab some data from an outside
server using `wget` on the command line.  

Make sure you are in the dc_workshop drectory first

    cd ~/dc_workshop
    wget http://reactomerelease.oicr.on.ca/download/archive/variant_calling.tar.gz

The file 'variant_calling.tar.gz' is what is commonly called a "tarball", which is 
a compressed archive similar to the .zip files we have seen before.  We can decompress 
this archive using the command below:

    tar xvf variant_calling.tar.gz

This will create a directory tree that contains some data (refence genome and fastq files).

<pre>
variant_calling
├── ref_genome
│   └── ecoli_rel606.fasta
├── run_variant_calling.sh
└── trimmed_fastq
    ├── SRR097977.fastq
    ├── SRR098026.fastq
    ├── SRR098027.fastq
    ├── SRR098028.fastq
    ├── SRR098281.fastq
    └── SRR098283.fastq
</pre>
  


