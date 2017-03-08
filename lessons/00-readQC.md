# Lesson QC of Sequence Read Data

Quality Control of NGS Data
===================

# Learning Objectives:
* Describe how the FastQ format encodes quality.
* Evaluate a FastQC report.
* Clean FastQ reads using Trimmommatic.
* Employ for loops to automate operations on multiple files.


## Details on the FASTQ format

Although it looks complicated  (and maybe it is), its easy to understand the [fastq](https://en.wikipedia.org/wiki/FASTQ_format) format with a little decoding. Some rules about the format include...

|Line|Description|
|----|-----------|
|1|Always begins with '@' and then information about the read|
|2|The actual DNA sequence|
|3|Always begins with a '+' and sometimes the same info in line 1|
|4|Has a string of characters which represent the quality scores; must have same number of characters as line 2|

so for example in our data set, one complete read is:
```
$ head -n4 SRR098281.fastq 
@SRR098281.1 HWUSI-EAS1599_1:2:1:0:318 length=35
CNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
+SRR098281.1 HWUSI-EAS1599_1:2:1:0:318 length=35
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
```
This is a pretty bad read. 

Notice that line 4 is:
```
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
```
As mentioned above, line 4 is a encoding of the quality. In this case, the code is the [ASCII](https://en.wikipedia.org/wiki/ASCII#ASCII_printable_code_chart) character table. According to the chart a '#' has the value 35 and '!' has the value 33 - **But these values are not actually the quality scores!** There are actually several historical differences in how Illumina and other players have encoded the scores. Here's the chart from wikipedia:

```
  SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS.....................................................
  ..........................XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX......................
  ...............................IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII......................
  .................................JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ......................
  LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL....................................................
  !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
  |                         |    |        |                              |                     |
 33                        59   64       73                            104                   126
  0........................26...31.......40                                
                           -5....0........9.............................40 
                                 0........9.............................40 
                                    3.....9.............................40 
  0.2......................26...31........41                              

 S - Sanger        Phred+33,  raw reads typically (0, 40)
 X - Solexa        Solexa+64, raw reads typically (-5, 40)
 I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
 J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
     with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator (bold) 
     (Note: See discussion above).
 L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)
 ```
 So using the Illumina 1.8 encoding, which is what you will mostly see from now on, our first C is called with a Phred score of 0 and our Ns are called with a score of 2. Read quality is assessed using the Phred Quality Score.  This score is logarithmically based and the score values can be interpreted as follows:

|Phred Quality Score |Probability of incorrect base call |Base call accuracy|
|:-------------------|:---------------------------------:|-----------------:|
|10	|1 in 10 |	90%|
|20	|1 in 100|	99%|
|30	|1 in 1000|	99.9%|
|40	|1 in 10,000|	99.99%|
|50	|1 in 100,000|	99.999%|
|60	|1 in 1,000,000|	99.9999%|

## FastQC
FastQC (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) provides a simple way to do some quality control checks on raw sequence data coming from high throughput sequencing pipelines. It provides a modular set of analyses which you can use to get a quick impression of whether your data has any problems of which you should be aware before doing any further analysis.

The main functions of FastQC are
* Import of data from BAM, SAM or FastQ files (any variant)
* Providing a quick overview to tell you in which areas there may be problems
* Summary graphs and tables to quickly assess your data
* Export of results to an HTML based permanent report
* Offline operation to allow automated generation of reports without running the interactive application


## Running FASTQC
### A. Stage your data

1. Create a working directory for your analysis
   
    ```bash
    $ cd
    # this command takes us to the home directory
    
    $ mkdir dc_workshop
    ```
2. Create three subdirectories

   ```bash
    mkdir dc_workshop/data
    mkdir dc_workshop/docs
    mkdir dc_workshop/results
```

  > The sample data we will be working with is in a hidden directory (placing a '.' in front of a directory name hides the directory. In the next step we will move some of those hidden files into our new dirctories to start our project. 
3. Move our sample data to our working (home) directory
   
   ```bash 
$ mv ~/.dc_sampledata_lite/untrimmed_fastq/ ~/dc_workshop/data/
```

###B. Run FastQC

1. Navigate to the initial fastq dataset
   
    ```bash
    $ cd ~/dc_workshop/data/untrimmed_fastq/
    ```
To run the fastqc program, we call it from its location in ``~/FastQC``.  fastqc will accept multiple file names as input, so we can use the *.fastq wildcard.
2. Run FastQC on all fastq files in the directory

    ```bash
    $ ~/FastQC/fastqc *.fastq
    ```
Now, let's create a home for our results
    ```bash
    $ mkdir ~/dc_workshop/results/fastqc_untrimmed_reads
    ```
3. Next, move the files there (recall, we are still in ``~/dc_workshop/data/untrimmed_fastq/``)
   ```bash 
    $ mv *.zip ~/dc_workshop/results/fastqc_untrimmed_reads/
    $ mv *.html ~/dc_workshop/results/fastqc_untrimmed_reads/
    ```


### C. Results

Let's examine the results in detail

1. Navigate to the results and view the directory contents

   ```bash
$ cd ~/dc_workshop/results/fastqc_untrimmed_reads/
$ ls
```
   
 > The zip files need to be unpacked with the 'unzip' program.  
2. Use unzip to unzip the FastQC results: 
   ```bash
$ unzip *.zip
```
Did it work? No, because 'unzip' expects to get only one zip file.  Welcome to the real world. We *could* do each file, one by one, but what if we have 500 files?  There is a smarter way. We can save time by using a simple shell 'for loop' to iterate through the list of files in *.zip. After you type the first line, you will get a special '>' prompt to type next next lines. You start with 'do', then enter your commands, then end with 'done' to execute the loop.
3. Build a ``for`` loop to unzip the files

   ```bash 
$ for zip in *.zip
> do
> unzip $zip
> done
```

  Note that, in the first line, we create a variable named 'zip'.  After that, we call that variable with the syntax $zip.  $zip is assigned the value of each item (file) in the list *.zip, once for each iteration of the loop.

This loop is basically a simple program.  When it runs, it will run unzip once for each file (whose name is stored in the $zip variable). The contents of each file will be unpacked into a separate directory by the unzip program.

The for loop is interpreted as a multipart command.  If you press the up arrow on your keyboard to recall the command, it will be shown like so:
   ```bash
    for zip in *.zip; do echo File $zip; unzip $zip; done
```

When you check your history later, it will help your remember what you did!

### D. Document your work

To save a record, let's cat all fastqc summary.txts into one full_report.txt and move this to ``~/dc_workshop/docs``. You can use wildcards in paths as well as file names.  Do you remember how we said 'cat' is really meant for concatenating text files?

```bash    
cat */summary.txt > ~/dc_workshop/docs/fastqc_summaries.txt
```


##How to clean reads using *Trimmomatic*
###A detailed explanation of features

Once we have an idea of the quality of our raw data, it is time to trim away adapters and filter out poor quality score reads. To accomplish this task we will use *Trimmomatic* [http://www.usadellab.org/cms/?page=trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic).

*Trimmomatic* is a java based program that can remove sequencer specific reads and nucleotides that fall below a certain threshold. *Trimmomatic* can be multithreaded to run quickly. 

Because *Trimmomatic* is java based, it is run using the command:

**_java -jar trimmomatic-0.32.jar_**

What follows this are the specific commands that tells the program exactly how you want it to operate. *Trimmomatic* has a variety of options and parameters:

* **_-threads_** How many processors do you want *Trimmomatic* to run with?
* **_SE_** or **_PE_** Single End or Paired End reads?
* **_-phred33_** or **_-phred64_** Which quality score do your reads have?
* **_SLIDINGWINDOW_** Perform sliding window trimming, cutting once the average quality within the window falls below a threshold.
* **_LEADING_** Cut bases off the start of a read, if below a threshold quality.
* **_TRAILING_** Cut bases off the end of a read, if below a threshold quality.
* **_CROP_** Cut the read to a specified length.
* **_HEADCROP_** Cut the specified number of bases from the start of the read.
* **_MINLEN_** Drop an entire read if it is below a specified length.
* **_TOPHRED33_** Convert quality scores to Phred-33.
* **_TOPHRED64_** Convert quality scores to Phred-64.

A generic command for *Trimmomatic* looks like this:

**java -jar trimmomatic-0.32.jar SE**

A complete command for *Trimmomatic* will look something like this:

**java -jar trimmomatic-0.32.jar SE -threads 4 -phred64 SRR_1056.fastq SRR_1056_trimmed.fastq ILLUMINACLIP:SRR_adapters.fa SLIDINGWINDOW:4:20**

This command tells *Trimmomatic* to run on a Single End file (``SRR_0156.fastq``, in this case), the output file will be called ``SRR_0156_trimmed.fastq``,  there is a file with Illumina adapters called ``SRR_adapters.fa``, and we are using a sliding window of size 4 that will remove bases with a phred score of below 20.


## Exercise - Running Trimmomatic

1. Go to the untrimmed fastq data location:

   ```bash
$ cd /home/dcuser/dc_workshop/data/untrimmed_fastq
```

The command line invocation for trimmomatic is more complicated.  This is where what you have been learning about accessing your command line history will start to become important.

The general form of the command is:

   ```bash
java -jar ~/Trimmomatic-0.32/trimmomatic-0.32.jar inputfile outputfile OPTION:VALUE...
```    
'java -jar' calls the Java program, which is needed to run trimmomatic, which lives in a 'jar' file (trimmomatic-0.32.jar), a special kind of java archive that is often used for programs written in the Java programing language.  If you see a new program that ends in '.jar', you will know it is a java program that is executed 'java -jar program name'.  The 'SE' argument is a keyword that specifies we are working with single-end reads.

The next two arguments are input file and output file names.  These are then followed by a series of options. The specifics of how options are passed to a program differ depending on the program. You will always have to read the manual of a new program to learn which way it expects its command-line arguments to be composed.


So, for the single fastq input file 'SRR098283.fastq', the command would be:
   ```bash
$ java -jar /home/dcuser/Trimmomatic-0.32/trimmomatic-0.32.jar SE SRR098283.fastq \
    SRR098283.fastq_trim.fastq SLIDINGWINDOW:4:20 MINLEN:20

    TrimmomaticSE: Started with arguments: SRR098283.fastq SRR098283.fastq_trim.fastq SLIDINGWINDOW:4:20 MINLEN:20
    Automatically using 2 threads
    Quality encoding detected as phred33
    Input Reads: 21564058 Surviving: 17030985 (78.98%) Dropped: 4533073 (21.02%)
    TrimmomaticSE: Completed successfully
```
So that worked and we have a new fastq file.

   ```bash
    $ ls SRR098283*
    SRR098283.fastq  SRR098283.fastq_trim.fastq
```

Now we know how to run trimmomatic but there is some good news and bad news.  One should always ask for the bad news first.  Trimmomatic only operates on one input file at a time and we have more than one input file.  The good news? We already know how to use a for loop to deal with this situation.

```bash
$ for infile in *.fastq
    >do
    >outfile=$infile\_trim.fastq
    >java -jar ~/Trimmomatic-0.32/trimmomatic-0.32.jar SE $infile $outfile SLIDINGWINDOW:4:20 MINLEN:20
    >done
```

Do you remember how the first word after "for" in the loop specifies a variable that is assigned the value of each item in the list in turn?  We can call it whatever we like.  This time it is called infile.  Note that the third line of this for loop is creating a second variable called outfile.  We assign it the value of $infile with '_trim.fastq' appended to it.  The '\' escape character is used so the shell knows that whatever follows \ is not part of the variable name $infile.  There are no spaces before or after the '='.







