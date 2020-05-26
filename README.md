# indoor-plants

## Script for chloroplast genome assembly

* This is a bash script I wrote to combine all the steps used in a chloroplast genome assembly. 
* The basic process is: get chloroplast-only reads out of all the reads, filter and downsample if required, assemble, polish with long reads, polish with short reads. 
* The script takes in long reads and short reads (R1 and R2). 
* You also need to specify a baits fasta file; that is - a file of sequences that are from a similar chloroplast genome to use to extract matching reads from other non-chloroplast reads (e.g. nuclear, mitochondria). I used three gene sequences from a similar species; the genes are likely to be fairly evenly spaced to allow matching reads to span the whole genome. 
* The script also needs a file of illumina adapter sequences for trimming. 
* The script activates a conda environment with the tools needed (this is bio.yml). A copy of this in in this repo. (It also makes a copy when you run it). 

Steps:
1/ extract chloroplast nanopore reads from all long reads
2/ keep a subset of these reads (longest, to x250 coverage) using filtlong
3/ assemble these with flye
4/ polish that assembly with the long reads themselves
5/ trim and filter illumina reads with fastp
6/ extract illumina chloroplast reads from all illumina reads
7/ downsample to X300 coverage with rasusa
8/ polish the flye assembly with the short reads, with pilon
9/ do an additional assembly with Unicycler, which uses both long and short reads
10/ do an additional assembly with miniasm, using long reads. 
11/ polish that assembly with minipolish and the long reads
12/ polish again with the short reads
13/ compare the assemblies with mummer - dnadiff, to get delta files (then view in assemblytics.com)
14/ calculate stats for all read sets and assemblies with seqkit stats. 







