---
layout: page
title: "Instructor Notes"
permalink: /guide/
---
# Instructor Notes for Wrangling Genomics

## Setting Up AWS

What to do during long command waits:

### File Downloads
The fastq files take about 15 minutes to download. This would be a good time to discuss the overall workflow of this lesson as illustrated by the graphic integrated on the page. It is recommended to start this lesson with the commmands to make and move to the /data/untrimmed-fastq directory and begin the download, and while files download, cover the "Bioinformatics Workflows" and "Starting with Data" texts. Beware that the last fastq file in the list takes the longest to download (~6-8 mins).

### Running FastQC
The FastQC analysis on all raw reads takes about 10 minutes to run. It is a good idea to have learners start this command and cover the FastQC background material and images while FastQC runs.

## Commands That Don't Run When Copied and Pasted
There are several commands that are example commands that will not run correctly if copy and pasted directly to the terminal. These commands serve as example commands and will need to be modified to fit each user. There is text around the commands outlining how they need to be changed, but it's helpful to be aware of them ahead of time as an instructor so you can set them up properly (if you can remember!).

#### scp Command to Download FastQC to local machines
In the FastQC section, learners will download FastQC output files in order to open '.html `.html` summary files on their local machines in a web browser. The scp command currently contains a public DNS (for example, `ec2-34-238-162-94.compute-1.amazonaws.com`), but this will need to be replaced with the public DNS of the machine used by each learner. The Public DNS for each learner will be the same one they use to log in. The password is again `data4Carp`. 

Command as is: 
~~~
scp dcuser@ec2-34-238-162-94.compute-1.amazonaws.com:~/dc_workshop/results/fastqc_untrimmed_reads/*.html ~/Desktop/fastqc_html
~~~

Command for learners to use:
~~~
scp dcuser@<Public_DNS>:~/dc_workshop/results/fastqc_untrimmed_reads/*.html ~/Desktop/fastqc_html
~~~

bwamem
trimmomatic
