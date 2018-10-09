---
title: "Variant Calling Workflow"
teaching: 35
exercises: 25
questions:
- "How do I find sequence variants between my sample and a reference genome?"
objectives:
- "Understand the steps involved in variant calling."
- "Describe the types of data formats encountered during variant calling."
- "Use command line tools to perform variant calling."
keypoints:
- "Bioinformatics command line tools are collections of commands that can be used to carry out bioinformatics analyses."
- "To use most powerful bioinformatics tools, you'll need to use the command line."
- "There are many different file formats for storing genomics data. It's important to understand what type of information is contained in each file, and how it was derived."
---

We mentioned before that we are working with files from a long-term evolution study of an *E. coli* population (designated Ara-3). Now that we have looked at our data to make sure that it is high quality, and removed low-quality base calls, we can perform variant calling to see how the population changed over time. We care how this population changed relative to the original population, *E. coli* strain REL606. Therefore, we will align each of our samples to the *E. coli* REL606 reference genome, and see what difference exist in our reads versus the genome.

# Alignment to a reference genome

![workflow_align](../img/variant_calling_workflow_align.png)

We perform read alignment or mapping to determine where in the genome our reads originated from. There are a number of tools to
choose from and, while there is no gold standard, there are some tools that are better suited for particular NGS analyses. We will be
using the [Burrows Wheeler Aligner (BWA)](http://bio-bwa.sourceforge.net/), which is a software package for mapping low-divergent
sequences against a large reference genome. 

The alignment process consists of two steps:

1. Indexing the reference genome
2. Aligning the reads to the reference genome


# Setting up

First we download the reference genome for *E. coli* REL606. Although we could copy or move the file with `cp` or `mv`, most genomics workflows begin with a download step, so we will practice that here. 

~~~
$ cd ~/dc_workshop
$ mkdir data/ref_genome
$ curl -L -o data/ref_genome/ecoli_rel606.fasta.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/985/GCA_000017985.1_ASM1798v1/GCA_000017985.1_ASM1798v1_genomic.fna.gz
$ gunzip data/ref_genome/ecoli_rel606.fasta.gz
~~~
{: .bash}

We will also download a set of trimmed FASTQ files to work with. These are small subsets of our real trimmed data, 
and will enable us to run our variant calling workflow quite quickly. 

~~~
# This is a place holder so the workflow will work. This will be updated with figshare links later
$ curl -L -o sub.tar.gz https://osf.io/ckf6w/download
$ tar xvf sub.tar.gz
$ mv sub/ ~/dc_workshop/data/trimmed_fastq_small
~~~
{: .bash}

You will also need to create directories for the results that will be generated as part of this workflow. We can do this in a single
line of code because `mkdir` can accept multiple new directory
names as input.

~~~
$ mkdir -p results/sam results/bam results/bcf results/vcf
~~~
{: .bash}


### Index the reference genome
Our first step is to index the reference genome for use by BWA. Indexing allows the aligner to quickly find potential alignment sites for query sequences in a genome, which saves time during alignment. Indexing the reference only has to be run once. The only reason you would want to create a new index is if you are working with a different reference genome or you are using a different tool for alignment.

~~~
$ bwa index data/ref_genome/ecoli_rel606.fasta
~~~
{: .bash}

While the index is created, you will see output something like this:

~~~
[bwa_index] Pack FASTA... 0.04 sec
[bwa_index] Construct BWT for the packed sequence...
[bwa_index] 1.05 seconds elapse.
[bwa_index] Update BWT... 0.03 sec
[bwa_index] Pack forward-only FASTA... 0.02 sec
[bwa_index] Construct SA from BWT and Occ... 0.57 sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa index data/ref_genome/ecoli_rel606.fasta
[main] Real time: 1.765 sec; CPU: 1.715 sec
~~~
{: .output}

### Align reads to reference genome

The alignment process consists of choosing an appropriate reference genome to map our reads against and then deciding on an 
aligner. We will use the BWA-MEM algorithm, which is the latest and is generally recommended for high-quality queries as it 
is faster and more accurate.

An example of what a `bwa` command looks like is below. This command will not run, as we do not have the files `ref_genome.fa`, `input_file_R1.fastq`, or `input_file_R2.fastq`.

~~~
$ bwa mem ref_genome.fasta input_file_R1.fastq input_file_R2.fastq > output.sam
~~~
{: .bash}

Have a look at the [bwa options page](http://bio-bwa.sourceforge.net/bwa.shtml). While we are running bwa with the default 
parameters here, your use case might require a change of parameters. *NOTE: Always read the manual page for any tool before using 
and make sure the options you use are appropriate for your data.*

We're going to start by aligning the reads from just one of the 
samples in our dataset (`SRRXXXXXXX.fastq`). Later, we'll be 
iterating this whole process on all of our sample files.

~~~
$ bwa mem data/ref_genome/ecoli_rel606.fasta data/trimmed_fastq_small/SRR2584863_1.trim.sub.fastq data/trimmed_fastq_small/SRR2584863_2.trim.sub.fastq > results/sam/SRR2584863.aligned.sam
~~~
{: .bash}

You will see output that starts like this: 

~~~
[M::bwa_idx_load_from_disk] read 0 ALT contigs
[M::process] read 78970 sequences (10000278 bp)...
[M::process] read 75986 sequences (10000157 bp)...
[M::mem_pestat] # candidate unique pairs for (FF, FR, RF, RR): (191, 37018, 12, 189)
[M::mem_pestat] analyzing insert size distribution for orientation FF...
[M::mem_pestat] (25, 50, 75) percentile: (428, 797, 1804)
[M::mem_pestat] low and high boundaries for computing mean and std.dev: (1, 4556)
[M::mem_pestat] mean and std.dev: (1087.05, 1032.12)
[M::mem_pestat] low and high boundaries for proper pairs: (1, 5932)
[M::mem_pestat] analyzing insert size distribution for orientation FR...
~~~
{: .output}


#### SAM/BAM format
The [SAM file](https://github.com/adamfreedman/knowyourdata-genomics/blob/gh-pages/lessons/01-know_your_data.md#aligned-reads-sam),
is a tab-delimited text file that contains information for each individual read and its alignment to the genome. While we do not 
have time to go in detail of the features of the SAM format, the paper by 
[Heng Li et al.](http://bioinformatics.oxfordjournals.org/content/25/16/2078.full) provides a lot more detail on the specification.

**The compressed binary version of SAM is called a BAM file.** We use this version to reduce size and to allow for *indexing*, which enables efficient random access of the data contained within the file.

The file begins with a **header**, which is optional. The header is used to describe source of data, reference sequence, method of
alignment, etc., this will change depending on the aligner being used. Following the header is the **alignment section**. Each line
that follows corresponds to alignment information for a single read. Each alignment line has **11 mandatory fields** for essential
mapping information and a variable number of other fields for aligner specific information. An example entry from a SAM file is 
displayed below with the different fields highlighted.

![sam_bam1](../img/sam_bam.png)


![sam_bam2](../img/sam_bam3.png)

We will convert the SAM file to BAM format using the `samtools` program with the `view` command and tell this command that the input is in SAM format (`-S`) and to output BAM format (`-b`): 

~~~
$ samtools view -S -b results/sam/SRR2584863.aligned.sam > results/bam/SRR2584863.aligned.bam
~~~
{: .bash}

~~~
[samopen] SAM header is present: 1 sequences.
~~~
{: .output}


### Sort BAM file by coordinates

Next we sort the BAM file using the `sort` command from `samtools`. `-o` tells the command where to write the output.

~~~
$ samtools sort -o results/bam/SRR2584863.aligned.sorted.bam results/bam/SRR2584863.aligned.bam 
~~~
{: .bash}

Our files are pretty small, so we won't see this output. If you run the workflow with larger files, you will see something like this:
~~~
[bam_sort_core] merging from 2 files...
~~~
{: .output}


SAM/BAM files can be sorted in multiple ways, e.g. by location of alignment on the chromosome, by read name, etc. It is important to be aware that different alignment tools will output differently sorted SAM/BAM, and different downstream tools require differently sorted alignment files as input.

You can use samtools to learn more about this bam file as well.

~~~
samtools flagstat results/bam/SRR2584863.aligned.sorted.bam
~~~
{: .bash}

This will give you the following statistics about your sorted bam file:

~~~
352019 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
2019 + 0 supplementary
0 + 0 duplicates
351867 + 0 mapped (99.96% : N/A)
350000 + 0 paired in sequencing
175000 + 0 read1
175000 + 0 read2
343112 + 0 properly paired (98.03% : N/A)
349706 + 0 with itself and mate mapped
142 + 0 singletons (0.04% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
~~~
{: .output}
## Variant calling

A variant call is a conclusion that there is a nucleotide difference vs. some reference at a given position in an individual genome
or transcriptome, often referred to as a Single Nucleotide Polymorphism (SNP). The call is usually accompanied by an estimate of 
variant frequency and some measure of confidence. Similar to other steps in this workflow, there are number of tools available for 
variant calling. In this workshop we will be using `bcftools`, but there are a few things we need to do before actually calling the 
variants.

![workflow](../img/variant_calling_workflow.png)

### Step 1: Calculate the read coverage of positions in the genome

Do the first pass on variant calling by counting read coverage with [bcftools](https://samtools.github.io/bcftools/bcftools.html). We will use the command `mpileup`. The flag `-O b` tells samtools to generate a bcf format output file, `-o` specifies where to write the output file, and `-f` flags the path to the reference genome:

~~~
$ bcftools mpileup -O b -o results/bcf/SRR2584863_raw.bcf \
-f data/ref_genome/ecoli_rel606.fasta results/bam/SRR2584863.aligned.sorted.bam 
~~~
{: .bash}

~~~
[mpileup] 1 samples in 1 input files
~~~
{: .output}

We have now generated a file with coverage information for every base.

### Step 2: Detect the single nucleotide polymorphisms (SNPs)

Identify SNPs using bcftools `call`. We have to specify ploidy with the flag `--ploidy`, which is one for the haploid *E. coli*. `-m` allows for multiallelic and rare-variant calling, `-v` tells the program to output variant sites only (not every site in the genome), and `-o` specifies where to write the output file:

~~~
$ bcftools call --ploidy 1 -m -v -o results/bcf/SRR2584863_variants.vcf results/bcf/SRR2584863_raw.bcf 
~~~
{: .bash}

### Step 3: Filter and report the SNP variants in variant calling format (VCF)

Filter the SNPs for the final output in VCF format, using `vcfutils.pl`:

~~~
$ vcfutils.pl varFilter results/bcf/SRR2584863_variants.vcf  > results/vcf/SRR2584863_final_variants.vcf
~~~
{: .bash}


## Explore the VCF format:

~~~
$ less -S results/vcf/SRR2584866_final_variants.vcf
~~~
{: .bash}

You will see the header (which describes the format), the time and date the file was
created, the version of bcftools that was used, the command line parameters used, and 
some additional information:

~~~
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##bcftoolsVersion=1.8+htslib-1.8
##bcftoolsCommand=mpileup -O b -o results/bcf/SRR2584863_raw.bcf -f data/ref_genome/ecoli_rel606.fasta results/bam/SRR2584863.aligned.sorted.bam
##reference=file://data/ref_genome/ecoli_rel606.fasta
##contig=<ID=CP000819.1,length=4629812>
##ALT=<ID=*,Description="Represents allele(s) other than observed.">
##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
##INFO=<ID=IDV,Number=1,Type=Integer,Description="Maximum number of reads supporting an indel">
##INFO=<ID=IMF,Number=1,Type=Float,Description="Maximum fraction of reads supporting an indel">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
##INFO=<ID=VDB,Number=1,Type=Float,Description="Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)",Version="3">
##INFO=<ID=RPB,Number=1,Type=Float,Description="Mann-Whitney U test of Read Position Bias (bigger is better)">
##INFO=<ID=MQB,Number=1,Type=Float,Description="Mann-Whitney U test of Mapping Quality Bias (bigger is better)">
##INFO=<ID=BQB,Number=1,Type=Float,Description="Mann-Whitney U test of Base Quality Bias (bigger is better)">
##INFO=<ID=MQSB,Number=1,Type=Float,Description="Mann-Whitney U test of Mapping Quality vs Strand Bias (bigger is better)">
##INFO=<ID=SGB,Number=1,Type=Float,Description="Segregation based metric.">
##INFO=<ID=MQ0F,Number=1,Type=Float,Description="Fraction of MQ0 reads (smaller is better)">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=ICB,Number=1,Type=Float,Description="Inbreeding Coefficient Binomial test (bigger is better)">
##INFO=<ID=HOB,Number=1,Type=Float,Description="Bias in the number of HOMs number (smaller is better)">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes for each ALT allele, in the same order as listed">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=DP4,Number=4,Type=Integer,Description="Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases">
##INFO=<ID=MQ,Number=1,Type=Integer,Description="Average mapping quality">
##bcftools_callVersion=1.8+htslib-1.8
##bcftools_callCommand=call --ploidy 1 -m -v -o results/bcf/SRR2584863_variants.bcf results/bcf/SRR2584863_raw.bcf; Date=Sat Jul  7 00:05:04 2018
~~~
{: .output}

Followed by information on each of the variations observed: 

~~~
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  results/bam/SRR2584863.aligned.sorted.bam
CP000819.1      9972    .       T       G       91      .       DP=4;VDB=0.0257451;SGB=-0.556411;MQ0F=0;AC=1;AN=1;DP4=0,0,0,4;MQ=60     GT:PL   1:121,0
CP000819.1      263235  .       G       T       85      .       DP=6;VDB=0.096133;SGB=-0.590765;RPB=1;MQB=1;BQB=1;MQ0F=0.166667;AC=1;AN=1;DP4=0,1,0,5;MQ=33     GT:PL   1:112,0
CP000819.1      281923  .       G       T       217     .       DP=10;VDB=0.774083;SGB=-0.662043;MQSB=0.974597;MQ0F=0;AC=1;AN=1;DP4=0,0,4,5;MQ=60       GT:PL   1:247,0
CP000819.1      433359  .       CTTTTTTT        CTTTTTTTT       64      .       INDEL;IDV=12;IMF=1;DP=12;VDB=0.477704;SGB=-0.676189;MQSB=1;MQ0F=0;AC=1;AN=1;DP4=0,1,3,8;MQ=60   G
CP000819.1      473901  .       CCGC    CCGCGC  228     .       INDEL;IDV=9;IMF=0.9;DP=10;VDB=0.659505;SGB=-0.662043;MQSB=0.916482;MQ0F=0;AC=1;AN=1;DP4=1,0,2,7;MQ=60   GT:PL   1
CP000819.1      648692  .       C       T       210     .       DP=10;VDB=0.268014;SGB=-0.670168;MQSB=0.916482;MQ0F=0;AC=1;AN=1;DP4=0,0,7,3;MQ=60       GT:PL   1:240,0
CP000819.1      1331794 .       C       A       178     .       DP=8;VDB=0.624078;SGB=-0.651104;MQSB=0.900802;MQ0F=0;AC=1;AN=1;DP4=0,0,3,5;MQ=60        GT:PL   1:208,0
CP000819.1      1733343 .       G       A       225     .       DP=11;VDB=0.992403;SGB=-0.670168;MQSB=1.00775;MQ0F=0;AC=1;AN=1;DP4=0,0,4,6;MQ=60        GT:PL   1:255,0
CP000819.1      2103887 .       ACAGCCAGCCAGCCAGCCAGCCAGCCAGCCAG        ACAGCCAGCCAGCCAGCCAGCCAGCCAGCCAGCCAGCCAGCCAGCCAGCCAGCCAG        56      .       INDEL;IDV=2;IMF=0.666667;
CP000819.1      2333538 .       AT      ATT     167     .       INDEL;IDV=7;IMF=1;DP=7;VDB=0.568173;SGB=-0.616816;MQSB=1.01283;MQ0F=0;
~~~
{: .output}

This is a lot of information, so let's take some time to make sure we understand our output.

The first few columns represent the information we have about a predicted variation. 

| column | info |
| ------- | ---------- |
| CHROM | contig location where the variation occurs | 
| POS | position within the contig where the variation occurs | 
| ID | a `.` until we add annotation information | 
| REF | reference genotype (forward strand) | 
| ALT | sample genotype (forward strand) | 
| QUAL | Phred-scaled probablity that the observed variant exists at this site (higher is better) |
| FILTER | a `.` if no quality filters have been applied, PASS if a filter is passed, or the name of the filters this variant failed | 

In an ideal world, the information in the `QUAL` column would be all we needed to filter out bad variant calls.
However, in reality we need to filter on multiple other metrics. 

The last two columns contain the genotypes and can be tricky to decode. 

| column | info |
| ------- | ---------- |
| FORMAT | lists in order the metrics presented in the final column | 
| results | lists the values associated with those metrics in order | 

For our file, the metrics presented are GT:PL:GQ. 

| metric | definition | 
| ------- | ---------- |
| GT | the genotype of this sample which for a diploid genome is encoded with a 0 for the REF allele, 1 for the first ALT allele, 2 for the second and so on. So 0/0 means homozygous reference, 0/1 is heterozygous, and 1/1 is homozygous for the alternate allele. For a diploid organism, the GT field indicates the two alleles carried by the sample, encoded by a 0 for the REF allele, 1 for the first ALT allele, 2 for the second ALT allele, etc. |
| PL | the likelihoods of the given genotypes |
| GQ | the Phred-scaled confidence for the genotype | 
| AD, DP | the depth per allele by sample and coverage |

The Broad Institute's [VCF guide](https://www.broadinstitute.org/gatk/guide/article?id=1268) is an excellent place
to learn more about VCF file format.

> ## Exercise
> 
> Use the `grep` and `wc` commands you've learned to assess how many variants are in the vcf file. 
>
>> ## Solution
>> 
>> ~~~
>> $ grep -v "#" results/vcf/SRR2584863_final_variants.vcf | wc -l
>> ~~~
>> {: .bash}
>> 
>> ~~~ 
>> 25
>> ~~~
>> {: .output}
>>
>> There are 25 variants in this file.
> {: .solution}
{: .challenge}

## Assess the alignment (visualization) - optional step

It is often instructive to look at your data in a genome browser. Visualisation will allow you to get a "feel" for 
the data, as well as detecting abnormalities and problems. Also, exploring the data in such a way may give you 
ideas for further analyses.  As such, visualization tools are useful for exploratory analysis. In this lesson we 
will describe two different tools for visualisation; a light-weight command-line based one and the Broad
Institute's Integrative Genomics Viewer (IGV) which requires
software installation and transfer of files.

In order for us to visualize the alignment files, we will need to index the BAM file using `samtools`:

~~~
$ samtools index results/bam/SRR2584863.aligned.sorted.bam
~~~
{: .bash}

### Viewing with `tview`

[Samtools](http://www.htslib.org/) implements a very simple text alignment viewer based on the GNU
`ncurses` library, called `tview`. This alignment viewer works with short indels and shows [MAQ](http://maq.sourceforge.net/) consensus. 
It uses different colors to display mapping quality or base quality, subjected to users' choice. Samtools viewer is known to work with an 130 GB alignment swiftly. Due to its text interface, displaying alignments over network is also very fast.

In order to visualize our mapped reads we use `tview`, giving it the sorted bam file and the reference file: 

~~~
$ samtools tview results/bam/SRR2584863.aligned.sorted.bam data/ref_genome/ecoli_rel606.fasta
~~~
{: .bash}

~~~
1         11        21        31        41        51        61        71        81        91        101       111       121
AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATAC
..................................................................................................................................
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, ..................N................. ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,........................
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, ..................N................. ,,,,,,,,,,,,,,,,,,,,,,,,,,,.............................
...................................,g,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,  ....................................   ................
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,....................................   ....................................      ,,,,,,,,,,
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,  ....................................  ,,a,,,,,,,,,,,,,,,,,,,,,,,,,,,,,     .......
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, .............................  ,,,,,,,,,,,,,,,,,g,,,,,    ,,,,,,,,,,,,,,,,,,,,,,,,,,,,
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,  ...........................T.......   ,,,,,,,,,,,,,,,,,,,,,,,c,          ......
......................... ................................   ,g,,,,,,,,,,,,,,,,,,,      ...........................
,,,,,,,,,,,,,,,,,,,,, ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, ,,,,,,,,,,,,,,,,,,,,,,,,,,,       ..........................
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,   ................................T..  ..............................   ,,,,,,
...........................       ,,,,,,g,,,,,,,,,,,,,,,,,   ....................................         ,,,,,,
,,,,,,,,,,,,,,,,,,,,,,,,,, ....................................  ...................................        ....
....................................  ........................  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,      ....
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,   ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
........................            .................................. .............................     ....
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,   ....................................        ..........................
...............................       ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, ....................................
...................................  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,  ..................................
.................................... ,,,,,,,,,,,,,,,,,,a,,,,,,,,,,,,,,,,,        ,,,,,,,,,,,,,,,,,,,,,,,,,
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,  ............................ ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
~~~
{: .output}

The first line of output shows the genome coordinates in our reference genome. The second line shows the reference
genome sequence. The third lines shows the consensus sequence determined from the sequence reads. A `.` indicates
a match to the reference sequence, so we can see that the consensus from our sample matches the reference in most
locations. That is good! If that wasn't the case, we should probably reconsider our choice of reference.

Below the horizontal line, we can see all of the reads in our sample aligned with the reference genome. Only 
positions where the called base differs from the reference are shown. You can use the arrow keys on your keyboard
to scroll or type `?` for a help menu. Type `Ctrl^C` to exit `tview`. 


### Viewing with IGV

[IGV](http://www.broadinstitute.org/igv/) is a stand-alone browser, which has the advantage of being installed locally and providing fast access. Web-based genome browsers, like [Ensembl](http://www.ensembl.org/index.html) or the [UCSC browser](https://genome.ucsc.edu/), are slower, but provide more functionality. They not only allow for more polished and flexible visualisation, but also provide easy access to a wealth of annotations and external data sources. This makes it straightforward to relate your data with information about repeat regions, known genes, epigenetic features or areas of cross-species conservation, to name just a few.

In order to use IGV, we will need to transfer some files to our local machine. We know how to do this with `scp`. 
Open a new tab in your terminal window and create a new folder. We'll put this folder on our Desktop for 
demonstration purposes, but in general you should avoide proliferating folders and files on your Desktop and 
instead organize files within a directory structure like we've been using in our `dc_workshop` directory.

~~~
$ mkdir ~/Desktop/files_for_igv
$ cd ~/Desktop/files_for_igv
~~~
{: .bash}

Now we will transfer our files to that new directory. Remember to replace the text between the `@` and the `:` 
with your AWS instance number. The commands to `scp` always go in the terminal window that is connected to your
local computer (not your AWS instance).

~~~
$ scp dcuser@ec2-34-203-203-131.compute-1.amazonaws.com:~/dc_workshop/results/bam/SRR2584863.aligned.sorted.bam ~/Desktop/files_for_igv
$ scp dcuser@ec2-34-203-203-131.compute-1.amazonaws.com:~/dc_workshop/results/bam/SRR2584863.aligned.sorted.bam.bai ~/Desktop/files_for_igv
$ scp dcuser@ec2-34-203-203-131.compute-1.amazonaws.com:~/dc_workshop/data/ref_genome/ecoli_rel606.fasta ~/Desktop/files_for_igv
$ scp dcuser@ec2-34-203-203-131.compute-1.amazonaws.com:~/dc_workshop/results/vcf/SRR2584863_final_variants.vcf ~/Desktop/files_for_igv
~~~
{: .bash}

You will need to type the password for your AWS instance each time you call `scp`. 

Next we need to open the IGV software. If you haven't done so already, you can download IGV from the [Broad Institute's software page](https://www.broadinstitute.org/software/igv/download), double-click the `.zip` file
to unzip it, and then drag the program into your Applications folder. 

1. Open IGV.
2. Load our reference genome file (`ecoli_rel606.fasta`) into IGV using the **"Load Genomes from File..."** option under the **"Genomes"** pull-down menu.
3. Load our BAM file (`SRR097977.aligned.sorted.bam`) using the **"Load from File..."** option under the **"File"** pull-down menu. 
4.  Do the same with our VCF file (`SRR097977_final_variants.vcf`).

Your IGV browser should look like the screenshot below:

![IGV](../img/igv-screenshot.png)

There should be two tracks: one coresponding to our BAM file and the other for our VCF file. 

In the **VCF track**, each bar across the top of the plot shows the allele fraction for a single locus. The second bar shows
the genotypes for each locus in each *sample*. We only have one sample called here so we only see a single line. Dark blue = 
heterozygous, Cyan = homozygous variant, Grey = reference.  Filtered entries are transparent.

Zoom in to inspect variants you see in your filtered VCF file to become more familiar with IGV. See how quality information 
corresponds to alignment information at those loci.
Use [this website](http://software.broadinstitute.org/software/igv/AlignmentData) and the links therein to understand how IGV colors the alignments.

Now that we've run through our workflow for a single sample, we want to repeat this workflow for our other five
samples. However, we don't want to type each of these individual steps again five more times. That would be very
time consuming and error-prone, and would become impossible as we gathered more and more samples. Luckily, we
already know the tools we need to use to automate this workflow and run it on as many files as we want using a
single line of code. Those tools are: wildcards, for loops, and bash scripts. We'll use all three in the next 
lesson. 

## Other Notes

#### Installing Software
 
It's worth noting here that all of the software we are using for
this workshop has been pre-installed on our remote computer. 
This saves us a lot of time - installing software can be a 
time-consuming and frustrating task - however, this does mean that
you won't be able to walk out the door and start doing these
analyses on your own computer. You'll need to install 
the software first. Look at the [setup instructions](http://www.datacarpentry.org/wrangling-genomics/setup/) for more information on installing these software packages.

#### BWA Alignment options
BWA consists of three algorithms: BWA-backtrack, BWA-SW and BWA-MEM. The first algorithm is designed for Illumina sequence 
reads up to 100bp, while the other two are for sequences ranging from 70bp to 1Mbp. BWA-MEM and BWA-SW share similar features such 
as long-read support and split alignment, but BWA-MEM, which is the latest, is generally recommended for high-quality queries as it 
is faster and more accurate. 
{: .callout}


> ## Multiple line commands
> 
> Some of the commands we ran in this lesson are burly! When typing a long 
> command into your terminal, you can use the `\` character
> to separate code chunks onto separate lines. This can make your code more readable.
>
{: .callout}
