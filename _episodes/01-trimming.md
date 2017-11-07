---
title: "Trimming and Filtering"
teaching: 0
exercises: 0
questions:
- "Is my data good enough?"
objectives:
- "Clean FASTQ reads using Trimmommatic."
- "Employ `for` loops to automate operations on multiple files."
keypoints:
- "First key point."
---

# How to clean reads using Trimmomatic

Some text here

## A. detailed explanation of features

Once we have an idea of the quality of our raw data, it is time to trim away adapters and filter out poor quality score reads. To accomplish this task we will use [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic).

Trimmomatic is a java based program that can remove sequencer specific reads and nucleotides that fall below a certain threshold. Trimmomatic can be multithreaded to run quickly. 

Because Trimmomatic is java based, it is run using the command:

~~~
$ java jar trimmomatic-0.32.jar
~~~
{: .bash}

What follows this are the specific commands that tells the program exactly how you want it to operate. Trimmomatic has a variety of options and parameters:

* `-threads` How many processors do you want *Trimmomatic* to run with?
* `SE` or `PE` Single End or Paired End reads?
* `-phred33` or `-phred64` Which quality score do your reads have?
* `SLIDINGWINDOW` Perform sliding window trimming, cutting once the average quality within the window falls below a threshold.
* `LEADING` Cut bases off the start of a read, if below a threshold quality.
* `TRAILING` Cut bases off the end of a read, if below a threshold quality.
* `CROP` Cut the read to a specified length.
* `HEADCROP` Cut the specified number of bases from the start of the read.
* `MINLEN` Drop an entire read if it is below a specified length.
* `TOPHRED33` Convert quality scores to Phred-33.
* `TOPHRED64` Convert quality scores to Phred-64.

A generic command for Trimmomatic looks like this:

~~~
$ java -jar trimmomatic-0.32.jar SE
~~~
{: .bash}

A complete command for Trimmomatic will look something like this:

~~~
java -jar trimmomatic-0.32.jar SE -threads 4 -phred64 SRR_1056.fastq SRR_1056_trimmed.fastq ILLUMINACLIP:SRR_adapters.fa SLIDINGWINDOW:4:20
~~~
{: .bash}

This command tells Trimmomatic to run on a Single End file (``SRR_0156.fastq``, in this case), the output file will be called ``SRR_0156_trimmed.fastq``,  there is a file with Illumina adapters called ``SRR_adapters.fa``, and we are using a sliding window of size 4 that will remove those bases if their phred score is below 20.

## Exercise - Running Trimmomatic

To begin, go to the untrimmed fastq data location:

~~~
$ cd /home/dcuser/dc_workshop/data/untrimmed_fastq
~~~
{: .bash}

The command line incantation for Trimmomatic is more complicated.  This is where what you have been learning about accessing your command line history will start to become important.

The general form of the command is:

~~~
$ java -jar ~/Trimmomatic-0.32/trimmomatic-0.32.jar inputfile outputfile OPTION:VALUE...
~~~
{: .bash}

`java -jar` calls the Java program, which is needed to run Trimmomatic, which lived in a "jar" file (trimmomatic-0.32.jar), a special kind of java archive that is often used for programs
written in the Java programing language.  If you see a new program that ends in `.jar`, you will know it is a java program that is executed `java -jar program name`.  The "SE" argument is a 
keyword that specifies we are working with single-end reads.

The next two arguments are input file and output file names.  These are then followed by a series of options. The specifics of how options are passed to a program are different depending on the 
program. You will always have to read the manual of a new program to learn which way it expects its command-line arguments to be composed.

So, for the single fastq input file `SRR098283.fastq`, the command would be:  

~~~
$ java -jar /home/dcuser/Trimmomatic-0.32/trimmomatic-0.32.jar SE SRR098283.fastq \
SRR098283.fastq_trim.fastq SLIDINGWINDOW:4:20 MINLEN:20
~~~
{: .bash}

~~~
TrimmomaticSE: Started with arguments: SRR098283.fastq SRR098283.fastq_trim.fastq SLIDINGWINDOW:4:20 MINLEN:20
Automatically using 2 threads
Quality encoding detected as phred33
Input Reads: 21564058 Surviving: 17030985 (78.98%) Dropped: 4533073 (21.02%)
TrimmomaticSE: Completed successfully
~~~
{: .output}

So that worked and we have a new fastq file.

~~~
$ ls SRR098283*
~~~
{: .bash}

~~~
SRR098283.fastq  SRR098283.fastq_trim.fastq
~~~
{: .output}

Now we know how to run Trimmomatic but there is some good news and bad news.  
One should always ask for the bad news first.  Trimmomatic only operates on 
one input file at a time and we have more than one input file.  The good news?
We already know how to use a `for` loop to deal with this situation.

~~~
$ for infile in *.fastq
> do
> outfile=$infile\_trim.fastq
> java -jar ~/Trimmomatic-0.32/trimmomatic-0.32.jar SE $infile $outfile SLIDINGWINDOW:4:20 MINLEN:20
> done
~~~
{: .bash}

Do you remember how the first specifies a variable that is assigned the value of each item in the list in turn?  We can call it whatever we like.  This time it is called `infile`. 
Note that the third line of this `for` loop is creating a second variable called `outfile`.  We assign it the value of `$infile` with `_trim.fastq` appended to it.
The `\` escape character is used so the shell knows that whatever follows `\` is not part of the variable name `$infile`.  There are no spaces before or after the `=`.
