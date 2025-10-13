#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=
#SBATCH --mem=
#SBATCH --time=
#SBATCH --array=0-99%10
#SBATCH --job-name=
#SBATCH --output=metaphlan_%A-%a.out

# set path
FILEPATH=`pwd -P`

# read list of sample names
readarray -t samples < $FILEPATH/sample_list.txt

# assign each sample to an array task
j=${samples[${SLURM_ARRAY_TASK_ID}]}

## -------------------------------------------------------------- ##
## FASTP - trim and filter reads
## -------------------------------------------------------------- ##

# make output directories
mkdir $FILEPATH/trimmed
mkdir $FILEPATH/read_qc

fastp \
--thread ${SLURM_CPUS_PER_TASK} \
-i $FILEPATH/data_raw/${j}_R1.fastq \
-I $FILEPATH/data_raw/${j}_R2.fastq \
-o $FILEPATH/trimmed/${j}_trimmed_R1.fastq \
-O $FILEPATH/trimmed/${j}_trimmed_R2.fastq \
--trim_poly_g \
--length_required 50 \
--qualified_quality_phred 25

mv fastp* $FILEPATH/read_qc

## -------------------------------------------------------------- ##
## BWA - remove human contamination
## -------------------------------------------------------------- ##

# make output directory
mkdir $FILEPATH/cleaned

bwa mem \
-t ${SLURM_CPUS_PER_TASK} \
$FILEPATH/data_reference/human_index/human_ref_seq \
$FILEPATH/trimmed/${j}_trimmed_R1.fastq $FILEPATH/trimmed/${j}_trimmed_R2.fastq \
| samtools sort -@16 -m16g - \
| samtools view -@16 -b -f 4 - \
| samtools fastq -@16 -1 $FILEPATH/cleaned/${j}_clean.1.fastq -2 $FILEPATH/cleaned/${j}_clean.2.fastq -s $FILEPATH/cleaned/${j}_singles.fastq -n

## -------------------------------------------------------------- ##
## METAPHLAN - identify taxonomic compositions
## -------------------------------------------------------------- ##

# make output directories
mkdir -p $FILEPATH/metaphlan_out/bowtie2_out

metaphlan \
$FILEPATH/cleaned/${j}_clean.1.fastq,$FILEPATH/cleaned/${j}_clean.2.fastq \
--bowtie2out $FILEPATH/metaphlan_out/bowtie2_out/${j}.bowtie2.bz2 \
--nproc ${SLURM_CPUS_PER_TASK} \
--input_type fastq \
-o $FILEPATH/metaphlan_out/${j}_profile.txt

## ------------------------------------------------------------- ##
## THE END
## ------------------------------------------------------------- ##
