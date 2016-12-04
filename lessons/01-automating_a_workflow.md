# Lesson

Shell scripts
===================

Learning Objectives:
-------------------
#### What's the goal for this lesson?

* Understand what a shell script is
* Learn how automate an analytical workflow


## What is a shell script?
A shell script is basically a text file that contains a list of commands
that are executed sequentially.  The commands in a shell script are the same
as you would use on the command line.

Once you have worked out the details and tested your commands in the shell, you can save them into a file so, the next time, you can automate the process with
a script.

The basic anatomy of a shell script is a file with a list of commands.
That is also the definition of pretty much any computer program.

```bash
#!/bin/bash

cd ~/dc_sample_data

for file in untrimmed_fastq/*.fastq
do
  echo "My file name is $file"
done
```

This looks a lot like the for loops we saw earlier.  In fact, it is no different, apart from using indentation and the lack of the '>' prompts; it's just saved in a text file. The line at the top ('#!/bin/bash') is commonly called the shebang line, which is a special kind of comment that tells the shell which program is to be used as the 'intepreter' that executes the code.  

In this case, the interpreter is bash, which is the shell environment we are working in. The same approach is also used for other scripting languages such as perl and python.  The shebang line is actually optionally unless you want to
make the script executable like a 'real' program.

## How to run a shell script
There are two ways to run a shell script. The first way is to specify the
interpreter (bash) and the name of the script.  By convention, shell script
use the .sh extension, though this is not enforced.

```bash
$ bash myscript.sh
My file name is untrimmed_fastq/SRR097977.fastq
My file name is untrimmed_fastq/SRR098026.fastq
```

The second way is a little more complicated to set up and requires the shebang line we talked about earlier.

The first step, which only needs to be done once, is to modify the 'permissions' of the text file so that the shell knows the file is executable.

```bash
$ chmod +x myscript.sh
```

After that, you can run the script as a regular program.

```bash
$ ./myscript.sh
$ bash myscript.sh 
My file name is untrimmed_fastq/SRR097977.fastq
My file name is untrimmed_fastq/SRR098026.fastq
```

The thing about running programs on the command line is that the shell may not know the location of your executables unless they are in the 'path' of known locations for programs.  So, you need to tell the shell the path to your script, which is './' if it is in the same directory.

****
**Exercise**
1) Use nano to save the code above to a script called myscript.sh
2) run the script
****


## A real shell script

Now, let's do something real.  First, recall the code from our our fastqc
workflow from this morning, with a few extra "echo" statements.

```bash
cd ~/dc_workshop/data/untrimmed_fastq/

echo "Running fastqc..."
~/FastQC/fastqc *.fastq
mkdir -p ~/dc_workshop/results/fastqc_untrimmed_reads

echo "saving..."
mv *.zip ~/dc_workshop/results/fastqc_untrimmed_reads/
mv *.html ~/dc_workshop/results/fastqc_untrimmed_reads/

cd ~/dc_workshop/results/fastqc_untrimmed_reads/

echo "Unzipping..."
for zip in *.zip
do
  unzip $zip
done

echo "saving..."
cat */summary.txt > ~/dc_workshop/docs/fastqc_summaries.txt
```


****
**Exercise**

1) Use nano to create a shell script using with the code above (you can copy/paste),
named read_qc.sh

2) Run the script

3) Bonus points: Use something you learned yesterday to save the output
of the script to a file while it is running.
****





