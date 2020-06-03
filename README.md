# indoor-plants

## Chloroplast genome assembly

* A bash script to assemble a chloroplast genome from long and short sequencing reads. 

## How to use

Example:

* Copy `plastid-assembly.sh` into the directory where you will run the analysis
* Activate your conda environment with the tools installed.  
* Run: `bash plastid-assembly.sh -a adapters.fasta -b baits.fasta R1.fq R2.fq nano.fq`

## Basic process

* Get chloroplast-only reads out of all the reads
* Filter and downsample if required
* Assemble
* Polish with long reads
* Polish with short reads

## Inputs

* Long reads, e.g. nanopore
* Short reads, e.g. Illumina R1 and R2
* An illumina adapters file, for trimming
* A baits.fasta file - sequences from a similar chloroplast genome to use to extract matching *chloroplast* reads from all the reads. I used three gene sequences from a similar species; the genes are likely to be fairly evenly spaced to allow matching reads to span the whole genome. 

## What the script does:

* extract chloroplast nanopore reads from all long reads
* keep a subset of these reads (longest, to x250 coverage) using filtlong
* assemble these with flye
* polish that assembly with the long reads themselves
* trim and filter illumina reads with fastp
* extract illumina chloroplast reads from all illumina reads
* downsample to X300 coverage with rasusa
* polish the flye assembly with the short reads, with pilon
* do an additional assembly with Unicycler, which uses both long and short reads
* do an additional assembly with miniasm, using long reads
* polish that assembly with minipolish and the long reads
* polish again with the short reads
* compare the assemblies with mummer - dnadiff, to get delta files (then view in http://www.assemblytics.com/)
* calculate stats for all read sets and assemblies with seqkit stats 

## Outputs

* Polished chloroplast genome assemblies using long and short reads, from Flye, Unicycler and Miniasm/Minipolish. 
* Upload the assembly fasta files to https://chlorobox.mpimp-golm.mpg.de/geseq.html if you want to annotate the genomes. 

## Tools

Installed with conda. 

```
minimap2
samtools
filtlong
flye
racon
fastp
rasusa
bwa
pilon
unicycler
miniasm
minipolish
mummer
seqkit
```

A conda env yml file is included this repo, with a record of the tools and versions used in the analysis, but this environment is platform-dependent (linux). 








