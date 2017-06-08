
# OSU Bioinformatics Workshop
***
## Genome Assembly Hands-on
Authors: Dana Brunson <dana.brunson@okstate.edu>, Peter Hoyt <peter.r.hoyt@okstate.edu> and Haibao Tang <tanghaibao@gmail.com>  

__Learning Outcomes:__ Get familiar with genome assemblers, pre-processing, reporting and validation. Exercise will be based on chromosomes of a mutant genotype of bakers’ yeast, as a practice of _de-novo_ genome assembly. Here we have the benefits of a reference genome to validate our assembly.

#### Data and directory structure ####

Log into the Cowboy Supercomputer (picture here):
https://hpcc.okstate.edu/content/logging-cowboy

Make a copy of all the data into your “/scratch/username” directory.
```
$ cd /scratch/username
$ cp /scratch/bioworkshop/mcbios.tar.gz .
$ tar xvf mcbios.tar.gz
$ cd mcbios/
$ ls

abyss data  results soap velvet 
```
Directory of folders you now have:
```
|-- abyss
|   `-- abyss31
|       `-- abyssk31.pbs
|-- data
|   |-- group1
|   |   |-- PE-350.1.fastq
|   |   |-- PE-350.2.fastq
|   |   `-- ref.fasta
|   |-- group2
|   |   |-- PE-350.1.fastq
|   |   |-- PE-350.2.fastq
|   |   `-- ref.fasta
|   |-- group3
|   |   |-- PE-350.1.fastq
|   |   |-- PE-350.2.fastq
|   |   `-- ref.fasta
|   |-- group4
|   |   |-- PE-350.1.fastq
|   |   |-- PE-350.2.fastq
|   |   `-- ref.fasta
|   `-- group5
|       |-- PE-350.1.fastq
|       |-- PE-350.2.fastq
|       `-- ref.fasta
|-- results
|   `-- quast.pbs
|-- soap
|   `-- soap31
|       |-- soap.config
|       `-- soapk31.pbs
`-- velvet
    `-- velvetk31.pbs
```


We will dive into each directory for each task:  fastqc, velvet, soap, abyss and nucmer in that order. Most folders contain a __submission script__ which includes the commands that we use for each task. It is always a good idea to use a script so you can modify parameters, and the script also serves as  a reminder to yourself later.
___
### Important notes before hands-on ###
Since we are using the Cowboy cluster, only very small tasks can be done directly on the login nodes.  For each longer activity, we will submit the jobs to the scheduler using “pbs scripts”.  These scripts include information for the scheduler as well as the commands to execute your job.
___
### Read processing ###
__Please Note: Those who have successfully completed “module2” of the Bioinformatics Applications/File Formats session (on Monday) can skip the Read Processing section.__


The data

```
    $ cd /scratch/USERNAME/mcbios/data/  OR cd data
    $ ls
    group1  group2  group3  group4  group5
```
__We will divide into 5 groups. Each group should work on a different chromosome in yeast.__

```$ cd group1
$ ls
ref.fasta  PE-350.1.fastq  PE-350.2.fastq
$ head PE-350.1.fastq
```
Our assembly will start using one Illumina library, PE-350, which is a paired end (fragment) library, divided into two datasets: .1.fastq and .2.fastq,  which are paired reads very close to 350 bp apart. Also, “ref.fasta” is the included reference sequence for comparison with your assembly.

#### Understand FASTQ format ####
__Learning objectives:__ 

What information is included? What does each line mean?
Four consecutive lines define ONE read in the fastq file

If we look at the first four lines of any .fastq file, we will see something like this:
```
@READ-ID
ATAGAGATAGAGAGAG (sequence read: in this example: 16 nucleotides)
+READ-ID (repeat the ID in first line, optional, can be empty)
;9;7;;.7;3933334 (quality identifiers/characters for each base; same length as sequence, in this example: 16 nucleotides)
```
So each fastq file has a total number of lines that is divisible by four. Can you use shell commands to count the number of sequence "reads" in a fastq file?
```
make the command line do all the arithmetic: 
$ expr $(cat PE-350.1.fastq |wc -l) / 4
```
Each individual base in a sequence read has quality identifier called PHRED score, typically between value 0-40 for Illumina.  The FASTQ file converts the numeric value to a character because they use less space (fewer bits). There are also two systems of such conversion, PHRED+33 and PHRED+64. PHRED+33 is used in new Illumina protocols. PHRED+64 is used also so be aware!

How does this work? 
“Phred quality scores  are defined as a property which is logarithmically related to the base-calling error probabilities .” -Wikipedia.

Or you can just use this chart: 

|Phred Quality Score |Probability of incorrect base call |Base call accuracy|
|:-------------------|:---------------------------------:|-----------------:|
|10	|1 in 10 |	90%|
|20	|1 in 100|	99%|
|30	|1 in 1000|	99.9%|
|40	|1 in 10,000|	99.99%|
|50	|1 in 100,000|	99.999%|
|60	|1 in 1,000,000|	99.9999%|

If your base has PHRED quality of 20 (i.e. 1% error), PHRED(20)+33=53, which corresponds to the character 5 using the reference chart below. Under PHRED+64 system, 20+64=84, which is T. When using the table below, can you guess if a sequence quality score is using PHRED+33, or PHRED+64 encoding? (hint: use the cheat sheet below which also shows a few other PHRED score schemes used)
```
  SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS.....................................................
  ..........................XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX......................
  ...............................IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII......................
  .................................JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ......................
  LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL....................................................
  !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
  |                         |    |        |                              |                     |
 33                        59   64       73                            104                   126
  0........................26...31.......40                                "Sanger"
                           -5....0........9.............................40 "Solexa"
                                 0........9.............................40 "Illumina 1.3+"
                                    3.....9.............................40 "Illumina 1.5+"
  0.2......................26...31........41                               "Illumina 1.8+"

 S - Sanger        Phred+33,  raw reads typically (0, 40)
 X - Solexa        Solexa+64, raw reads typically (-5, 40)
 I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
 J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
 L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)
 ```
Go back to our FASTQ files,  find out answers to these questions:

1. Library PE-350 contains a total of ______ reads. (hint: count # of lines using the command line by typing:  `wc -l PE-350.1.fastq` and divide by 4). Leave out the quotes. 
2. The base quality in my data is encoded in ________ (Phred+33 or Phred+64 ?).

3. How good is your sequencing run? Use FASTQC
#### Running FASTQC ####
Fastqc runs very quickly, so we can run this on the login nodes.  
`$ module load fastqc`

For info on fastqc (read this later):
`$ fastqc -h`

Now perform the fastqc quality control
`$ fastqc PE-350.1.fastq`

(__protip:__ do this again to the PE-350.2.fastq file by using the up arrow
to bring up the last command, then use arrow keys to move over to change 
1 to a 2 and press enter.)

`$ fastqc PE-350.2.fastq`

fastqc puts the results in the same folder as the data (in this case data/group1/)

__Instructions to download your fastqc results are at this link.__ 

FASTQC generates a HTML report (which means you can open it in any web browser) for each FASTQ file you run.  Use the WINSCP program to find the folder in /scratch/username/mcbios/data/group __X__ which holds the results of your analysis.  

Use [WinSCP](https://winscp.net/eng/download.php) to transfer the PE-350.1_fastqc folder to your desktop and double-click on  file __fastqc_report.html__ (on your desktop) to open it in a browser. For more information and to see examples of what bad data look like, read the FASTQC manual when you have time. 


#### Report on Library PE-350 fastQC stats: #### 
Read count_____, Read length ___ 
Duplication level _____, Adapters? _____ (hint: look for“over-represented sequences” in report)

Answer the same questions for PE-350.2.fastq

### ASSEMBLY ###

__Learning Objectives:__ Not all assemblers are the same or give the same results.

Which assembler should I use?

It is (in general) a good idea to run several assemblers and compare. We are using VELVET, SOAPdenovo and ABYSS in this study. All of these assemblers are [de-bruijn](https://en.wikipedia.org/wiki/De_Bruijn_graph) (K-mer) graph based assemblers.

#### VELVET ####
```
$ cd ../../velvet/
```
Understand the script file: `velvetk31.pbs` - We will use the text-file editor “nano” to examine this file, which submits our data to the assemblers. 
```
$ nano -w velvetk31.pbs
```
__Using nano:__  To exit press Ctrl-x (hold down the ctrl key and press x.)  If you have made changes, it will ask if you want to save the changes you made.  Press 'y' to for yes to save. Then just press 'Enter' to accept the filename `velvetk31.pbs`. 

Be sure to change the line:    
"export GROUPNUMBER=1"    
to __your__ group number.  __You will need to do this for all the submission scripts today.__

The first section of our `velvetk31.pbs` script sets up the scheduler requests (the queueing system), and environment variables. The  assembler commands are `velveth` and `velvetg`. 

What do these commands do?  See the [velvet website here.](https://www.ebi.ac.uk/~zerbino/velvet/)  

If you understand the script and have changed the groupnumber, press `ctrl-x`, then save the file and __submit it to the queue__ with the ‘qsub’ command as follows.
```
$ qsub velvetk31.pbs 
```
What was printed on the screen after you pressed enter?
 
`<jobidnumber>.mgmt1` 

The log file, `velvetk31.pbs.o<jobidnumbmer>` is very useful. 

Options to examine the jobidnumber file:
```
$ tail velvetk31.pbs.o<jobidnumber>  (HINT: Try using TAB-completion to enter your commands!)
OR
$ less velvetk31.pbs.o<jobidnumber> 
(q to quit)
```
#### Challenges: ####
Find these information from the log file:

+ Paired end library insert size: _____________   
+ Standard deviation ______________
+ contig stats: n50 ________  max _______ total __________ reads used ____________ (did it use all the reads?)

+ Where are the contigs stored? Save the results!
```
$ cp velvet31/contigs.fa ../results/velvet31.fasta
```
+ Does the assembly get better  if I use a different K-mer size? (you can try this later)

To try different kmers, first copy your pbs script:
```
$ cp velvetk31.pbs velvetk21.pbs
$ nano -w velvetk21.pbs
```
and change K to a different value, currently 31. Try 21, 25.
then submit this new job:
```
$ qsub velvetk21.pbs
```
To see kmer coverage differences in a histogram format, we can use the following Perl script on the velvet stats file in each of your results folders: 
```
$ module load velvet
$ velvet-estimate-exp_cov.pl velvet31/stats.txt
    10 |      1 | **********
<snip>
Predicted expected coverage: 18
velvetg parameters: -exp_cov 18 -cov_cutoff 0
```
CONGRATULATIONS on your first assembly!
___

#### SOAPdenovo2  [website](http://soap.genomics.org.cn/soapdenovo.html) ####

Now we will compare our data using the SOAPdenovo assembler. 

First change to the soap31 subdirectory
```
$ cd ../soap/soap31
$ ls 
soap.config  soapk31.pbs
```
SOAPdenovo is different as it uses a configuration file “soap.config”, where we tell SOAP what our data are. We include detailed info below so we can move ahead with the computation step.







`[LIB]`
Calls the soap denovo command library into action

`avg_ins=350`
This value indicates the average insert size of this library or the peak value position in the insert size distribution figure.

`reverse_seq=0`
This option takes value 0 or 1. It tells the assembler if the read sequences need to be complementarily reversed. Illumima GA produces two types of paired-end libraries: a) forward-reverse, generated from fragmented DNA ends with typical insert size less than 500 bp; b) reverse-forward, generated from circularizing libraries with typical insert size greater than 2 Kb. The parameter "reverse_seq" should be set to indicate this: 0, forward-reverse; 1, reverse-forward.

`asm_flags=3`
This indicator decides in which part(s) the reads are used. It takes value 1 (only contig assembly), 2 (only scaffold assembly), 3 (both contig and scaffold assembly), or 4 (only gap closure).

`rank=1`
It takes integer values and decides in which order the reads are used for scaffold assembly. Libraries with the same "rank" are used at the same time during scaffold assembly.
 
Challenge: How do you set library rank?

SOAPdenovo will use the paired-end libraries with insert sizes from smaller to larger to construct scaffolds. Libraries with the ___same rank___ would be used at the same time. For example, in a dataset of a human genome, we set five ranks for five libraries with insert size 200-bp, 500-bp, 2-Kb, 5-Kb and 10-Kb, separately. It is desired that the pairs in each rank provide adequate physical coverage of the genome.

`q1=../../data/group1/PE-350.1.fastq`
Pair-end file 1 (with quality scores)

`q2=../../data/group1/PE-350.2.fastq`
Pair-end file 2 (with quality scores)

First edit the soap.config file to the correct data for your group:

```$ nano -w soap.config
Change these two lines from group1 to your group:
q1=../../data/group1/PE-350.1.fastq
q2=../../data/group1/PE-350.2.fastq
```
Understand `soapk31.pbs`, it is important to know SOAPdenovo has 4 steps: __pregraph__, __contig__, __map__, and __scaff__ are ‘step 1’ and correspond to making a K-mer graph, contigging, mapping the reads back, and scaffolding are steps 2-4 respectively.  These steps can be run separately or all together as we are doing here.  See the command line options section of the manual for more information. 

At this point you don’t need to change the submit script, if you want to look at it:
```
$ more soapk31.pbs
```
and to submit the soapdenovo assembly:
```
$ qsub soapk31.pbs
```
What printed on the screen? 
ans: <jobid>.mgmt1

The __‘soapk31.pbs.o<jobid>’__ is very useful. 

Challenge: Use the ‘more’ command to find these information from the log file:
+ Paired end library insert size: _____________   
+ Standard deviation ______________
+ contig stats: n50 ________  max length contig _______ total __________ 

Where are the results stored?
FASTA file: soap31.scafSeq

Now SOAP is done. Save the results!
```
$ cp soap31.scafSeq ../../results/soap31.fasta
```
Does the assembly get better  if I use a different K-mer size?

To run a different kmer e.g ‘21’, create a new soap directory and copy the submit and config files into it. Notice that we are changing the names of the files while copying them to reflect the new kmer we are testing:
```
$ cd ..
$ pwd (make sure you’re in the soap directory)
$ mkdir soap21
$ cp soap31/soapk31.pbs soap21/soapk21.pbs
$ cp soap31/soap.config soap21/.      (that dot at the end is necessary)
$ cd soap21
$ nano -w soapk21.pbs 
```
Now modify `soapk21.pbs` and change K  from 31 to a different value. Try 21,  then 25.  After submitting your job:
```
qsub soapk21.pbs
```
Be sure to copy the result to your appropriate results folder.
```
$ cp soap21.scafSeq ../../results/soap21.fasta
```
Using what you’ve learned, ___do the 25 kmer value___. 
___

#### ABYSS ####
```
$ cd ../../abyss/abyss31
```

We change to a subdirectory because abyss puts all it’s output into the current working directory.

Understand `abyssk31.pbs`. ABYSS is simple as it has just one command to run the entire pipeline. First edit the `abyssk31.pbs` to change GROUPNUMBER to __your__ group number.
```
$ nano -w abyssk31.pbs
$ qsub abyssk31.pbs
```
##### K-mer coverage #####

There is one file of interest here: `coverage.hist` as generated by ABYSS. This is the K-mer histogram. It is a two-column file with 1st column (coverage), 2nd column (# of K-mers). 

Plot this K-mer histogram however you want (EXCEL, R, etc.). Only the first 50 lines are very useful, so we will help you shrink the file to include just the first 50 lines:

To look at the first 50 lines, enter:
```
$ head -n 50 coverage.hist
```
To create a new file that only includes these first 50 lines. :
```
$ head -50 coverage.hist >short.hist
```

EMAIL the ‘short.hist’ file to yourself to look at it in a spreadsheet.
```
$ mail -a short.hist -r youremail@wherever.com  youremail@wherever.com
(enter subject, enter, ctrl-d to send)
```
To determine the final assembled size of the genome, use `grep` at the command line to extract that information from the output file:
```
$ grep ^Assembled abyssk31.pbs.o<jobidnumber>
Assembled 563522 k-mer in 245 contigs.
```
for more information: https://groups.google.com/forum/#!topic/abyss-users/RdR6alqM7e8

#### Challenges ####
To determine the highest K-mer coverage in the histogram use:
```
$ grep median abyssk31.pbs.o<jobidnumber>
```
and put your answer in this box ___________

You can estimate genome size just based on this K-mer histogram (why do you want to do that?)

__Genome_size = Total_Kmers / peak_Kmer_coverage__

(hint: To calculate Total_Kmers, Create a third column by multiplying column 1 and column2, and then sum the numbers of the 3rd column)

Use the formula to estimate the genome size base on K-mer histogram: ______________ bp

The `abyssk31.pbs.o<jobid>` is very useful. Find these information from the log file:
contig stats: n50 ________  max _______ total __________ (compare to the estimate based on K-mer)
scaffold stats: n50 ________  max _______ total __________ 

Where are the results stored?

Now ABYSS is done. Save the results and rename it!
```
$ cp abyss31-scaffolds.fa ../../results/abyss31.fasta
```
Does the assembly get better  if I use a different K-mer size? 

__Changing the K-mer option for assemblies:__  Luckily assemblers run fast, run two different additional K-mer options using Abyss. Why? Because the current value K=31 may not be the best! Try 21, and 25 for most assemblers, they need odd numbers, so 24 won’t work.  

First create a new directory:
```
$ cd ..
```
(you should be in the ‘abyss directory now, to check use “print working directory”)
```
$ pwd  
/scratch/username/mcbios/abyss
$ mkdir abyss21
$ cp abyss31/abyssk31.pbs abyss21/abyssk21.pbs
$ cd abyss21
$ nano -w abyssk21.pbs
```
 change K to a different value, currently 31. Try 21, 25 and submit:
```
$ qsub abyssk21.pbs
```
Whan it’s done, your output files will automatically have the new K-mer in their names. Copy the result to your results folder
```
$ cp abyss21-scaffolds.fa ../../results/abyss21.fasta
```

__Reporting:__

A _contig_ is a contiguous length of genomic sequence. A _scaffold_ is composed of ordered contigs and gaps. By far the most widely used statistics for describing the quality of a genome assembly are its scaffold and contig N50s. 

A contig N50 is calculated by first ordering every contig by length from longest to shortest. Next, starting from the longest contig, the lengths of each contig are summed, until this running sum equals one-half of the total length of all contigs in the assembly. The contig N50 of the assembly is the length of the shortest contig in this list. 

(Insert image here)

The scaffold N50 is calculated in the same fashion but uses scaffolds rather than contigs. The longer the scaffold N50 is, the better the assembly is. However, it is important to keep in mind that a poor assembly that has forced unrelated reads and contigs into scaffolds can have an erroneously large N50.

N50 statistic is a metric of the length of a set of sequences. N50 is the contig length such that using equal or longer contigs produces half the bases (http://en.wikipedia.org/wiki/N50_statistic).

```
$ cd ../../results
$ ls
```
This directory contains the results from all the programs. If you have a different assembly (e.g different K-mers) using the same assembler, name them differently, for example: `abyss25.fasta`,  `abyss31.fasta`,  `soap31.fasta`,  `velvet31.fasta`

__Note:__ Abyss has a command to get the scaffold statistic (N50) of an assembly.

Use can use the command `abyss-fac`, which scans the contigs lengths and outputs the N50 for your assembly. We will cover K-mer comparisons later but thought you might like to know:
```
$ module load abyss
$ abyss-fac velvet31.fasta
```
___
### Validation ###
Make sure that you validate the results before releasing it. Some assemblies may appear to have large contigs and scaffolds, but are wrong. ___Check the assembly!___
#### QUAST ###
One software package for evaluating assemblies is quast.py [website](http://quast.sourceforge.net/quast)
Make sure you are in `results` directory (otherwise, `cd results`)
```
$ nano -w quast.pbs
```
The file looks like this:
```
#!/bin/bash
#
#PBS -q express
#PBS -j oe
#PBS -l nodes=1:ppn=12
#PBS -l walltime=1:00:00
cd $PBS_O_WORKDIR

module load quast

export GROUPNUMBER=1
export DATADIR=../data/group${GROUPNUMBER}

quast.py --gene-finding  `ls *.fasta`  -o quastresults  -R ${DATADIR}/ref.fasta
```
Remember to change your GROUPNUMBER.  It will analyze all the *.fasta files in your results directory

submit quast.pbs
```
$ qsub quast.pbs
```
when finished look at the output file to check for errors
```
$ less quast.pbs.o<jobid>  (protip: TAB autocomplete so you don’t have to type in the jobid)
```
press “q” to exit 

options:  zip the quast directory and mail the whole thing to yourself
```
$ zip -r quast.zip quastresults
$ mail -a quast.zip -r <youremailaddress> <youremailaddress>
(remember to hit ctrl-d to send)
```
Extract the zip file and just double-click on the `report.html` file.  Note that this html (web page) is *interactive*.

Other useful output from QUAST:

`report.pdf`: contains similar analyses and plots as the html file (double-click to open)

`alignment.svg`: contains the contig alignment plot (double-click to open)

+ Blue blocks: correctly aligned. Boundaries agree (within 2 kbp on each side, contigs are larger than 10 kbp) in at least half of the assemblies
+ Green blocks: correctly aligned. Boundaries don’t agree. 
+ Orange blocks: misassembled. Boundaries agree in at least half of the assemblies.
+ Red blocks: misassembled. Boundaries don’t agree. 

#### Challenge ####
Which assembly has better length statistics? 
Check out the [Quast manual](http://quast.bioinf.spbau.ru/manual.html) to get answers.

Best assembly: _________________

Worst assembly: _________________

How does your assembly compare structurally to the reference? 

Align contigs to your best assembly, change the red filename below. 
```
$ cd ../
$ mkdir nucmer
$ cd nucmer
$ module load bio_apps
$ nucmer ../results/abyss31.fasta ../data/group1/ref.fasta
```
Visualize it in dot plot:
```
$ mummerplot out.delta --postscript --layout
$ ps2pdf out.ps
$ mail -a out.pdf -r <youremailaddress> <youremailaddress>
(remember to hit ctrl-d to send)
```
Repeat the above procedure on the worst assembly.  How different is it ? ____________________


Background: about the data used in the exercise:

The data used in this exercise comes from a mutant yeast, using a novel method to generate the mutant. Our data comes from an individual called `MUTATOR4`.

Exercise originally presented at OSU by Dr. Haibao Tang, JCVI. Modified with extraordinary effort by Dr. Dana Brunson, OSU HPCC, and presented by Dr. Peter R. Hoyt

The paper:
http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3364565/

The data (we partitioned the reads to chromosomes so that assembly ran faster in workshop):
http://www.ncbi.nlm.nih.gov/sra/DRX001304


This concludes the exercise. We’ll continue SNP analysis with this dataset in another module.


EXTRAS:
http://kmergenie.bx.psu.edu/  “kmergenie”
See some slides here: http://ged.msu.edu/angus/tutorials-2013/files/2013-june-18-msu.pdf
Another de novo assembly tutorial: http://www.cbs.dtu.dk/courses/27626/Exercises/denovo_exercise.php
(uses quake & jellyfish)
good outline: http://en.wikibooks.org/wiki/Next_Generation_Sequencing_(NGS)/De_novo_assembly
github repostiory with lots of slides etc: https://github.com/lexnederbragt/denovo-assembly-tutorial


