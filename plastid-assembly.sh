#!/usr/bin/env bash

#A script to assemble a chloroplast genome

#............................................................................
# How to use

# activate your conda environment with the tools needed
# bash plastid-assembly.sh -a adapters -b baits R1 R2 nano

#............................................................................
# Defaults

set -e #exit if a command exits with non zero status
script=$(basename $0) #script name less file path
adapters=""
baits=""
genome_size=160000
threads=16

#............................................................................
# Functions

function msg {
  echo -e "$*"
}
# $* is args passed in, -e means interpret backslashes

function err {
  echo "Error: $*" 1>&2
}
# redirects stdout to stderr

function banner {
  printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' -
}

function msg_banner {
  printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' -
  msg "$*"
  printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' -
}

function usage {
  banner
  msg "$script\n   Assemble a chloroplast genome with long and short reads."
  msg "Author:\n  Anna Syme, anna.syme@gmail.com"
  msg "Usage:\n   $script [options] R1 R2 nano"
  msg "Parameters:"
  msg "   Illumina R1 reads      R1.fq.gz"
  msg "   Illumina R2 reads      R2.fq.gz"
  msg "   Nanopore reads         nano.fq.gz"
  msg "Options:"
  msg "   -h                     Show this help"
  msg "   -t NUM                 Number of threads (default=16)"
  msg "   -g NUM                 genome size in bp (default=160000)"
  msg "   -a FILE                illumina adapter sequences file name"
  msg "   -b FILE                bait sequences file name"
  msg "Example:"
  msg "   $script -t 8 R1.fq.gz R2.fq.gz minion.fq.gz"
  banner
  exit 1
  #exits script
}
#...........................................................................
# Parse the command line options

#loops through, sets variable or runs function
#have to add in a colon after the flag, if it's looking for an arg
#e.g. t: means flag t is set, then look for the arg (eg 32) to set $threads to
# : at the start disables default error handling for this bit
#instead it will label an odd flag as a ?
#the extra : at the end means it will label missing args as :

while getopts ':ht:g:a:b::' opt ; do
  case $opt in
    h)
      usage
      ;;
    t)
      threads=$OPTARG
      ;;
    g)
      genome_size=$OPTARG
      ;;
    a)
      adapters=$OPTARG
      ;;
    b)
      baits=$OPTARG
      ;;
    \?)
      echo "Invalid option '=$OPTARG'"
      exit 1
      ;;
    :)
      echo "Option '-$OPTARG' requires an argument"
      exit 1
      ;;
  esac
done

shift $((OPTIND-1))
#remove all options that has been parsed by getopts
#so that $1 refers to next arg passed to script
#e.g. a positional arg

if [ $# -ne 3 ]; then
  msg "\n **Please provide three input parameters** \n"
  usage
fi
#if number of pos params is not 3, print usage msg

#...........................................................................
#start
msg "\n"
msg_banner "now running $script"

#...........................................................................
#check inputs

#give variable names to the inputs
R1_raw=$1
R2_raw=$2
nano_raw=$3

msg "This script will use:"
msg "   Illumina reads R1:   $R1_raw"
msg "   Illumina reads R2:   $R2_raw"
msg "   Nanpore reads:       $nano_raw"
msg "   Bait sequences:      $baits"
msg "   Adapter sequences:   $adapters"
msg "   Genome size:         $genome_size"
msg "   Threads:             $threads"

#...........................................................................

conda env export --name plastidenv > plastidenv.yml
#saves conda env with tools and versions

#...........................................................................
#these read files will be created during run
R1_fastp=R1_fastp.fq.gz
R2_fastp=R2_fastp.fq.gz
R1_cp=R1_cp.fq.gz
R2_cp=R2_cp.fq.gz
R1_cp_subset=R1_cp_subset.fq.gz
R2_cp_subset=R2_cp_subset.fq.gz
nano_cp=nano_cp.fq.gz
nano_cp_long=nano_cp_long.fq.gz

#...........................................................................
#these assemblies will be created during run
assembly_flye=assembly_flye.fasta
assembly_flye_racon1=assembly_flye_racon1.fasta
assembly_flye_racon2=assembly_flye_racon2.fasta
assembly_flye_racon_pilon1=assembly_flye_racon_pilon1.fasta
assembly_unicycler=assembly_unicycler.fasta
assembly_miniasm=assembly_miniasm.fasta
assembly_miniasm_minipolished=assembly_miniasm_minipolished.fasta
assembly_miniasm_minipolished_pilon1=assembly_miniasm_minipolished_pilon1.fasta

#...........................................................................
msg_banner "now extracting chloroplast nanopore reads from all reads"

minimap2 -a -x map-ont -t $threads $baits $nano_raw | 
samtools fastq -0 $nano_cp -n -F 4 -
#map the reads to a reference set of chloroplast genes (baits)
#baits chosen are rbcL, matK and ndhF from Acacia ligulata
#this pipes the output to samtools fastq,
#which extracts the fastq reads from the alignment
#the flag -F 4 means exclude unmapped reads (i.e., non cp reads)

#...........................................................................
msg_banner "now keeping only the longest nanopore cp reads"

filtlong --length_weight 1 --mean_q_weight 0 --window_q_weight 0 \
--target_bases 40000000 $nano_cp | gzip > $nano_cp_long

##keeps set of longest reads to required cov (X200) based on genome size (160000)
##or cov (X250) = 40000000
##read score based on length only due to relative score weights
##alternative: filter by length and quality e.g.
##filtlong --target_bases 32000000 reads_nano_cp.fq.gz \
##| gzip > reads_nano_cp_filtered.fq.gz

#...........................................................................
msg_banner "now assembling nanopore cp reads"

flye --nano-raw $nano_cp_long --genome-size $genome_size \
--out-dir flye-out --threads $threads
#using uncorrected reads as read correction can lead to errors
#other flye options had little effect here - keep haplotpyes, meta, trestle
cp flye-out/assembly.fasta $assembly_flye

#...........................................................................
msg_banner "now polishing flye assembly with long reads"

#round 1. make overlaps file, then run racon
minimap2 -x map-ont -t $threads $assembly_flye $nano_cp_long \
| gzip > overlaps1.paf.gz
racon --threads $threads $nano_cp_long overlaps1.paf.gz \
$assembly_flye > $assembly_flye_racon1

#round 2
minimap2 -x map-ont -t $threads $assembly_flye_racon1 $nano_cp_long \
| gzip > overlaps2.paf.gz
racon --threads $threads $nano_cp_long overlaps2.paf.gz \
$assembly_flye_racon1 > $assembly_flye_racon2
#here, further rounds of polishing made little difference
#option to add in medaka polishing here

#...........................................................................
msg_banner "now trimming and filtering raw illumina reads"

fastp --in1 $R1_raw --out1 $R1_fastp --in2 $R2_raw --out2 $R2_fastp --verbose \
--adapter_fasta $adapters --n_base_limit 3 --length_required 130 \
--average_qual 35 --thread $threads

#filter illumina reads for quality and trim adapters
#--n_base_limit 3 - discard read/pair if > 3Ns
#--length_required 130
#--average_qual 35 - want very high qual for polishing

#...........................................................................
msg_banner "now extracting illumina chloroplast reads from all reads"

minimap2 -a -x sr $assembly_flye_racon2 $R1_fastp $R2_fastp \
| samtools fastq -1 $R1_cp -2 $R2_cp -F 0x4 -f 0x2 -
#extract cp reads only by mapping to long-read assembly
#chose assembly made with longest nanopore reads
#needs end dash for stdin
#more info https://broadinstitute.github.io/picard/explain-flags.html
#samtools flags: -F4 exclude unmapped reads, -f2 include properly paired reads

#...........................................................................
msg_banner "now creating subset of illumina chloroplast reads"

rasusa -i $R1_cp -i $R2_cp --coverage 300 \
--genome-size $genome_size -o $R1_cp_subset -o $R2_cp_subset
#downsample to required coverage, x300

#...........................................................................
msg_banner "now polishing flye assembly with illumina"

#round 1pilon polish
bwa index $assembly_flye_racon2
bwa mem -t $threads $assembly_flye_racon2 $R1_cp_subset $R2_cp_subset \
| samtools sort > flye_aln1.bam
samtools index flye_aln1.bam
samtools faidx $assembly_flye_racon2
pilon --genome $assembly_flye_racon2 --frags flye_aln1.bam \
--output assembly_flye_racon_pilon1 \
--fix bases --mindepth 0.5 --changes --threads $threads --verbose
#fix bases, not contig breaks in case that makes incorrect breaks

#round 2 pilon polish
#made no difference but left in as option
#bwa index $assembly_flye_racon_pilon1
#bwa mem -t $threads $assembly_flye_racon_pilon1 $R1_cp_subset $R2_cp_subset \
#| samtools sort > flye_aln2.bam
#samtools index flye_aln2.bam
#samtools faidx $assembly_flye_racon_pilon1
#pilon --genome $assembly_flye_racon_pilon1 --frags flye_aln2.bam \
#--output assembly_flye_racon_pilon2 \
#--fix bases --mindepth 0.5 --changes --threads $threads --verbose

#...........................................................................
msg_banner "now running unicycler assembler"

unicycler -1 $R1_cp_subset -2 $R2_cp_subset -l $nano_cp_long -o unicycler \
--threads $threads --no_rotate --keep 2
#using fastp filtered illumina reads and longest nano cp reads
#--keep 2 will keep final files but also SAM
#--no_rotate means don't rotate replicons to certain start pos
cp unicycler/assembly.fasta $assembly_unicycler

#...........................................................................
msg_banner "now running miniasm assembler"

#map reads to themselves - miniasm needs this as input
minimap2 -x ava-ont $nano_cp_long $nano_cp_long | gzip -1 > reads_overlaps.paf.gz
#assemble with miniasm - needs reads and overlaps (the paf file)
miniasm -f $nano_cp_long reads_overlaps.paf.gz > miniasm.gfa
#convert gfa to fasta file of unitigs
awk '/^S/{print ">"$2"\n"$3}' miniasm.gfa > $assembly_miniasm

#...........................................................................
msg_banner "now polishing miniasm assembly with long reads"

#minipolish, uses racon and attempts to circularise
minipolish -t $threads $nano_cp_long miniasm.gfa > minipolished.gfa
awk '/^S/{print ">"$2"\n"$3}' minipolished.gfa > $assembly_miniasm_minipolished

#alternative: use racon polish only
#racon polish round 1
#minimap2 -x map-ont $assembly_miniasm $nano_cp_long | gzip > overlaps1.paf.gz
#racon --threads $threads $nano_cp_long overlaps1.paf.gz \
#$assembly_miniasm > $assembly_miniasm_racon1

#racon polish round 2
#minimap2 -x map-ont $assembly_miniasm_racon1 $nano_cp_long | gzip > overlaps2.paf.gz
#racon --threads $threads $nano_cp_long overlaps2.paf.gz \
#$assembly_miniasm_racon1 > $assembly_miniasm_racon2

#...........................................................................
msg_banner "now polishing miniasm assembly with short reads"

#pilon polish round 1
bwa index $assembly_miniasm_minipolished
bwa mem -t $threads $assembly_miniasm_minipolished $R1_cp_subset $R2_cp_subset \
| samtools sort > mini_pilon_aln1.bam
samtools index mini_pilon_aln1.bam
samtools faidx $assembly_miniasm_minipolished

pilon --genome $assembly_miniasm_minipolished --frags mini_pilon_aln1.bam \
--output assembly_miniasm_minipolished_pilon1 \
--fix bases --mindepth 0.5 --changes --threads $threads --verbose
#fix bases, not contig breaks in case that makes incorrect breaks

#pilon polish round 2
#option. made no difference here
#bwa index $assembly_miniasm_racon_pilon1
#bwa mem -t $threads $assembly_miniasm_racon_pilon1 $R1_cp_subset $R2_cp_subset \
#| samtools sort > mini_pilon_aln2.bam
#samtools index mini_pilon_aln2.bam
#samtools faidx $assembly_miniasm_racon_pilon1
#pilon --genome $assembly_miniasm_racon_pilon1 --frags mini_pilon_aln2.bam \
#--output assembly_miniasm_racon_pilon2 \
#--fix bases --mindepth 0.5 --changes --threads $threads --verbose

#...........................................................................
msg_banner "now comparing assemblies with dnadiff"

#compare assemblies with mummer dnadiff
#$assembly_unicycler
#$assembly_flye_racon_pilon1
#$assembly_miniasm_minipolished_pilon1
#also include the uncorrected, plain miniasm unitigs
#$assembly_miniasm

dnadiff -p miniasm-unicycler $assembly_miniasm $assembly_unicycler
dnadiff -p miniasm-flye $assembly_miniasm $assembly_flye_racon_pilon1
dnadiff -p miniasm-miniasm_corrected \
$assembly_miniasm $assembly_miniasm_minipolished_pilon1
#see file.report for info
#then view delta file in assemblytics.com

#...........................................................................
msg_banner "now calculating stats for reads and assemblies"

seqkit stats $nano_raw $nano_cp $nano_cp_long \
-Ta > nano_read_stats.tsv

seqkit stats $R1_raw $R2_raw $R1_fastp $R2_fastp $R1_cp $R2_cp \
$R1_cp_subset $R2_cp_subset \
-Ta > illumina_read_stats.tsv

seqkit stats \
$assembly_flye \
$assembly_flye_racon1 $assembly_flye_racon2 $assembly_flye_racon_pilon1 \
$assembly_unicycler \
$assembly_miniasm \
$assembly_miniasm_minipolished $assembly_miniasm_minipolished_pilon1 \
-Ta > assembly_stats.tsv

#option: canu assembler
#canu -p canu -d canu genomeSize=160000 -nanopore $reads

#...........................................................................
msg_banner "Script finished!"
