---
title: "Trimming and Filtering"
teaching: 30
exercises: 25
questions:
- "How can I get rid of sequence data that doesn't meet my quality standards?"
objectives:
- "Clean FASTQ reads using Trimmommatic."
- "Select and set multiple options for command-line bioinformatics tools."
- "Write `for` loops with two variables."
keypoints:
- "The options you set for the command-line tools you use are important!"
- "Data cleaning is an essential step in a genomics workflow."
---

# Cleaning Reads

In the previous episode, we took a high-level look at the quality
of each of our samples using FastQC. We vizualized per-base quality
graphs showing the distribution of read quality at each base across
all reads in a sample and extracted information about which samples
fail which quality checks. We know that all of our samples failed at
least one of the quality metrics used by FastQC. This doesn't mean,
though, that our samples should be thrown out! It's very common to
have some reads within a sample,
or some positions (near the begining or end of reads) across all
reads that are low
quality and should be discarded. We will use a program called
[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) to
filter poor quality reads and trim poor quality bases from our samples.

## Trimmomatic Options

Trimmomatic has a variety of options to trim your reads. If we run the commands, we can see some of our options.

~~~
$ trimmomatic
~~~
{: .bash}

Which will give you the following output:
~~~
Usage: 
       PE [-version] [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-summary <statsSummaryFile>] [-quiet] [-validatePairs] [-basein <inputBase> | <inputFile1> <inputFile2>] [-baseout <outputBase> | <outputFile1P> <outputFile1U> <outputFile2P> <outputFile2U>] <trimmer1>...
   or: 
       SE [-version] [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-summary <statsSummaryFile>] [-quiet] <inputFile> <outputFile> <trimmer1>...
   or: 
       -version
~~~
{: .output}

This output shows us that we must first specify whether we have paired end (`PE`) or single end (`SE`) reads.
Next, we specify what flag we would like to run. For example, you can specify `threads` to indicate the number
processors on your computer that you want Trimmomatic to use. These flags are not necessary, but they can give
you more control over the command. The flags are followed by positional arguments, meaning the order in which you 
specify them is important. In paired end mode, Trimmomatic expects the two input files, and then the names of the
output files. These files are described below.

| option    | meaning |
| ------- | ---------- |
|  <inputFile1>  | Input reads to be trimmed. Typically the file name will contain an `_1` or `_R1` in the name.|
| <inputFile2> | Input reads to be trimmed. Typically the file name will contain an `_2` or `_R2` in the name.|
|  <outputFile1P> | Output file that contains surviving pairs from the `_1` file. |
|  <outputFile1U> | Output file that contains orphaned reads from the `_1` file. |
|  <outputFile2P> | Output file that contains surviving pairs from the `_2` file.|
|  <outputFile2U> | Output file that contains orphaned reads from the `_2` file.|

The last thing trimmomatic expects to see is the trimming parameters:

| step   | meaning |
| ------- | ---------- |
| `ILLUMINACLIP` | Perform adapter removal |
| `SLIDINGWINDOW` | Perform sliding window trimming, cutting once the average quality within the window falls below a threshold. |
| `LEADING`  | Cut bases off the start of a read, if below a threshold quality.  |
|  `TRAILING` |  Cut bases off the end of a read, if below a threshold quality. |
| `CROP`  |  Cut the read to a specified length. |
|  `HEADCROP` |  Cut the specified number of bases from the start of the read. |
| `MINLEN`  |  Drop an entire read if it is below a specified length. |
|  `TOPHRED33` | Convert quality scores to Phred-33.  |
|  `TOPHRED64` |  Convert quality scores to Phred-64. |

We will use only a few of these options and trimming steps in our
analysis. It is important to understand the steps you are using to
clean your data. For more information about the Trimmomatic arguments
and options, see [the Trimmomatic manual](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf).


However, a complete command for Trimmomatic will look something like this:

~~~
$ trimmomatic PE -threads 4 SRR_1056_1.fastq SRR_1056_2.fastq  \
              SRR_1056_1.trimmed.fastq SRR_1056_1un.trimmed.fastq \
              SRR_1056_2.trimmed.fastq SRR_1056_2un.trimmed.fastq
              ILLUMINACLIP:SRR_adapters.fa SLIDINGWINDOW:4:20
~~~
{: .bash}

In this example, we've told Trimmomatic:

| code   | meaning |
| ------- | ---------- |
| `PE` | that it will be taking a single end file as input |
| `-threads 4` | to use four computing threads to run (this will spead up our run) |
| `SRR_1056_1.fastq` | the first input file name |
| `SRR_1056_2.fastq` | the second input file name |
| `SRR_1056_1.trimmed.fastq` | the output file for surviving pairs from the `_1` file |
| `SRR_1056_1un.trimmed.fastq` | the output file for orphaned reads from the `_1` file |
| `SRR_1056_2.trimmed.fastq` | the output file for surviving pairs from the `_2` file |
| `SRR_1056_2un.trimmed.fastq` | the output file for orphaned reads from the `_2` file |
| `ILLUMINACLIP:SRR_adapters.fa`| to clip the Illumina adapters from the input file using the adapter sequences listed in `SRR_adapters.fa` |
|`SLIDINGWINDOW:4:20` | to use a sliding window of size 4 that will remove bases if their phred score is below 20 |

## Running Trimmomatic

Now we will run Trimmomatic on our data. To begin, navigate to your `untrimmed_fastq` data directory:

~~~
$ cd ~/dc_workshop/data/untrimmed_fastq
~~~
{: .bash}

We are going to run Trimmomatic on one of our paired-end samples. We
will use a sliding window of size 4 that will remove bases if their
phred score is below 20 (like in our example above). We will also
discard any reads that do not have at least 20 bases remaining after
this trimming step.

~~~
$ trimmomatic PE SRR098283.fastq SRR098283.fastq_trim.fastq SLIDINGWINDOW:4:20 MINLEN:25
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

> ## Exercise
>
> Use the output from your Trimmomatic command to answer the
> following questions.
>
> 1) What percent of reads did we discard from our sample?
> 2) What percent of reads did we keep?
>
>> ## Solution
>> 1) 21.02%
>> 2) 78.98%
> {: .solution}
{: .challenge}

You may have noticed that Trimmomatic automatically detected the
quality encoding of our sample. It is always a good idea to
double-check this or to enter the quality encoding manually.

We can confirm that we have our output file:

~~~
$ ls SRR098283*
~~~
{: .bash}

~~~
SRR098283.fastq  SRR098283.fastq_trim.fastq
~~~
{: .output}

The output file is also a FASTQ file. It should be smaller than our
input file because we've removed reads. We can confirm this:

~~~
$ ls SRR098283* -l -h
~~~
{: .bash}

~~~
-rw-r--r-- 1 dcuser dcuser 3.9G Jul 30  2015 SRR098283.fastq
-rw-rw-r-- 1 dcuser dcuser 3.0G Nov  7 23:10 SRR098283.fastq_trim.fastq
~~~
{: .output}


We've just successfully run Trimmomatic on one of our FASTQ files!
However, there is some bad news. Trimmomatic can only operate on
one sample at a time and we have more than one sample. The good news
is that we can use a `for` loop to iterate through our sample files
quickly!

~~~
$ for infile in *.fastq
> do
> outfile="${infile}"_trim.fastq
> trimmomatic PE "${infile}" "${outfile}" SLIDINGWINDOW:4:20 MINLEN:20
> done
~~~
{: .bash}

The new part in our `for` loop is the line:

~~~
> outfile="${infile}"_trim.fastq
~~~
{: .bash}

`infile` is the first variable in our loop and takes the value
of each of the FASTQ files in our directory. `outfile` is the
second variable in our loop and is defined by adding `_trim.fastq` to
the end of our input file name. Use `{}` to wrap the variable so that `_trim.fastq` will
not be interpreted as part of the variable name. In addition, quoting the shell variables is
a good practice AND necessary if your variables have spaces in them. For more, check [Bash Pitfalls](http://mywiki.wooledge.org/BashPitfalls).
There are no spaces before or after the `=`.

Go ahead and run the for loop. It should take a few minutes for
Trimmomatic to run for each of our six input files. Once it's done
running, take a look at your directory contents.

~~~
$ ls
~~~
{: .bash}

~~~
SRR097977.fastq		    SRR098027.fastq_trim.fastq	SRR098283.fastq
SRR097977.fastq_trim.fastq  SRR098028.fastq		SRR098283.fastq_trim.fastq
SRR098026.fastq		    SRR098028.fastq_trim.fastq	SRR098283.fastq_trim.fastq_trim.fastq
SRR098026.fastq_trim.fastq  SRR098281.fastq
SRR098027.fastq		    SRR098281.fastq_trim.fastq
~~~
{: .output}

If you look very closely, you'll see that you have three files for the
`SRR098283` sample. This is because we already had the `SRR098283.fastq_trim.fastq` file in our directory when we started
our `for` loop (because we had run Trimmomatic on just that one file already).
Our `for` loop included this file in our list of `.fastq` files and
created a new output file named `SRR098283.fastq_trim.fastq_trim.fastq`, which is the result of
running Trimmomatic on our already trimmed file. `SRR098283.fastq_trim.fastq` and `SRR098283.fastq_trim.fastq_trim.fastq` should be identical. If you look at your Trimmomatic output in the terminal window, you will see:

~~~
TrimmomaticSE: Started with arguments: SRR098283.fastq_trim.fastq SRR098283.fastq_trim.fastq_trim.fastq SLIDINGWINDOW:4:20 MINLEN:20
Automatically using 2 threads
Quality encoding detected as phred33
Input Reads: 17030985 Surviving: 17030985 (100.00%) Dropped: 0 (0.00%)
TrimmomaticSE: Completed successfully
~~~
{: .output}

This shows that when we re-trimmed our trimmed file, no new reads were
dropped. This is a good thing!

> ## Exercise
> Earlier we looked at the first read in our `SRR098026.fastq` file and
> saw that it was very poor quality.
>
> ~~~
> $ head -n4 SRR098026.fastq
> ~~~
> {: .bash}
>
> ~~~
> @SRR098026.1 HWUSI-EAS1599_1:2:1:0:968 length=35
> NNNNNNNNNNNNNNNNCNNNNNNNNNNNNNNNNNN
> +SRR098026.1 HWUSI-EAS1599_1:2:1:0:968 length=35
> !!!!!!!!!!!!!!!!#!!!!!!!!!!!!!!!!!!
> ~~~
> {: .output}
>
> After filtering out bad reads, what is the first remaining read for
> this sample? What does its quality look like?
>
>> ## Solution
>> ~~~
>> $ head -n4 SRR098026.fastq_trim.fastq
>> ~~~
>> {: .bash}
>>
>> ~~~
>> @SRR098026.342 HWUSI-EAS1599_1:2:1:3:655 length=35
>> GGATNGGCCTTGTATTTATGATTCTCNGAGTCTGT
>> +SRR098026.342 HWUSI-EAS1599_1:2:1:3:655 length=35
>> BB@B!B@AACBBABCCCCBBBBBB@@!B?B<ABB@
>> ~~~
>> {: .output}
>>
>> Comparing this with our quality scale:
>>
>> ~~~
>> Quality encoding: !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI
>>                   |         |         |         |         |
>> Quality score:    0........10........20........30........40
>> ~~~
>> {: .output}
>>
>> We can see that the scores are mostly in the 30+ range. This is
>> pretty good.
> {: .solution}
{: .challenge}

We've now completed the trimming and filtering steps of our quality
control process! Before we move on, let's move our trimmed FASTQ files
to a new subdirectory within our `data/` directory. We can also remove
our extra, double-trimmed file for the `SRR098283` sample.

~~~
$ cd ~/dc_workshop/data/untrimmed_fastq
$ mkdir ../trimmed_fastq
$ rm SRR098283.fastq_trim.fastq_trim.fastq
$ mv *_trim* ../trimmed_fastq
$ cd ../trimmed_fastq
$ ls
~~~
{: .bash}

~~~
SRR097977.fastq_trim.fastq  SRR098028.fastq_trim.fastq
SRR098026.fastq_trim.fastq  SRR098281.fastq_trim.fastq
SRR098027.fastq_trim.fastq  SRR098283.fastq_trim.fastq
~~~
{: .output}

> ## Bonus Exercise (Advanced)
>
> Now that we've quality controled our samples, they should perform
> better on the quality tests run by FastQC. Go ahead and re-run
> FastQC on your trimmed FASTQ files and visualize the HTML files
> to see whether your per base sequence quality is higher after
> trimming.
>
>> ## Solution
>>
>> In your AWS terminal window do:
>>
>> ~~~
>> $ ~/FastQC/fastqc ~/dc_workshop/data/trimmed_fastq
>> ~~~
>> {: .bash}
>>
>> In a new tab in your terminal do:
>>
>> ~~~
>> $ mkdir ~/Desktop/fastqc_html/trimmed
>> $ scp dcuser@ec2-34-203-203-131.compute-1.amazonaws.com:~/dc_workshop/data/trimmed_fastq/*.html ~/Desktop/fastqc_html/trimmed
>> $ open ~/Desktop/fastqc_html/trimmed/*.html
>> ~~~
>> {: .bash}
>>
>> Remember to replace everything between the `@` and `:` in your scp
>> command with your AWS instance number.
>>
>> Before trimming, one of the sequences gave a warning and another
>> failed the per base sequence quality test. After filtering, all
>> sequences pass that test.
> {: .solution}
{: .challenge}
