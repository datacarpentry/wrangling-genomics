# Metadata table notes

## Blount et al. 2012

Genomic analysis of a key innovation in an experimental Escherichia coli population
http://dx.doi.org/10.1038/nature11514
supplementary table 1: "historical Ara-3 clones subjected to whole genome sequencing"
Notes:
+ changed clade cit+ to C3+ or C3+H to match notation in Leon et al. 2018
+ used information in supplementary table 1: "historical Ara-3 clones subjected to whole genome sequencing"

## Tenaillon et al. 2016

Tempo and mode of genome evolution in a 50,000-generation experiment
http://dx.doi.org/10.1038/nature18959
supplementary data 1: https://media.nature.com/original/nature-assets/nature/journal/v536/n7615/extref/nature18959-s1.xlsx

## Leon et al. 2018

Innovation in an E. coli evolution experiment is contingent on maintaining adaptive potential until competition subsides
https://doi.org/10.1371/journal.pgen.1007348
S1 Table. Genome sequencing of E. coli isolates from the LTEE population.
Clade designations describe placement in the phylogenetic tree of all sequenced strains from the population and relative to key evolutionary transitions in this population: UC, Unsuccessful Clade; C1, Clade 1; C2, Clade 2; C3, Clade 3; C3+, Clade 3 Cit+; C3+H, Clade 3 Cit+ hypermutator.
https://doi.org/10.1371/journal.pgen.1007348.s006

https://tracykteal.github.io/introduction-genomics/01-intro-to-dataset.html states that genome sizes are not real data -- I haven't added these in yet. 
Ecoli_metadata.csv downloaded from http://www.datacarpentry.org/R-genomics/data/Ecoli_metadata.csv 

### Other changes

+ make all column headers lower case 
+ There was conflicting information for three strains. I chose to represent Blount et al. 2012 in the master sheet:
	+ ZDB99 is recorded as C2 in Leon et al. 2018, but as C1 in Blount et al. 2012. 
	+ ZDB30 is recorded as C3+ (cit+) by Leon et al. 2018, but as C3 (cit-) in Blount et al. 2012
	+ ZDB143 is recorded in Leon et al. 2018 as C2, but as Cit+ in Blount et al. 2012
+ When data is missing, I kept the cell blank