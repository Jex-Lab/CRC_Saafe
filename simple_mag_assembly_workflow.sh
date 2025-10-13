#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=
#SBATCH --mem=
#SBATCH --time=
#SBATCH --array=0-99%10
#SBATCH --job-name=
#SBATCH --output=mags_%A-%a.out

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
## METASPADES - assemble reads into contigs
## -------------------------------------------------------------- ##

# make output directory
mkdir $FILEPATH/assembly

metaspades.py \
-1 $FILEPATH/cleaned/${j}_clean.1.fastq.gz \
-2 $FILEPATH/cleaned/${j}_clean.2.fastq.gz \
-o $FILEPATH/assembly/${j} \
--threads ${SLURM_CPUS_PER_TASK}

## -------------------------------------------------------------- ##
## METAQUAST - assembly QC
## -------------------------------------------------------------- ##

# make output directory
mkdir $FILEPATH/asembly_stats

metaquast \
-o $FILEPATH/assembly_stats/${j} \
--fast \
--threads ${SLURM_CPUS_PER_TASK} \
$FILEPATH/assembly/${j}/contigs.fasta

## -------------------------------------------------------------- ##
## METABINNER - create binning input files
## -------------------------------------------------------------- ##

# make output directory
mkdir $FILEPATH/binning

bash /path_to_script/gen_coverage_file.sh \
-a $FILEPATH/assembly/${j}/contigs.fasta \
-o $FILEPATH/binning/${j} \
-t ${SLURM_CPUS_PER_TASK} \
-m ${SLURM_MEM_PER_NODE} \
-l 1000 \
-f .1.fastq.gz \
-r .2.fastq.gz \
$FILEPATH/trimmed/${j}_clean*

python /path_to_script/gen_kmer.py \
$FILEPATH/assembly/${j}/contigs.fasta \
1000 4

## -------------------------------------------------------------- ##
## METABINNER - bin contigs into mags
## -------------------------------------------------------------- ##

bash run_metabinner.sh \
-a $FILEPATH/assembly/${j}/contigs.fasta \
-o $FILEPATH/binning/${j} \
-d $FILEPATH/binning/${j}/coverage_profile_f1k.tsv \
-k $FILEPATH/assembly/${j}/contigs_kmer_4_f1000.csv \
-p /path_to_script \
-t ${SLURM_CPUS_PER_TASK}

## -------------------------------------------------------------- ##
## DAStool - write mags to fasta (also performs dereplication)
## -------------------------------------------------------------- ##

DAS_Tool \
--bins $FILEPATH/binning/${j}/metabinner_res/metabinner_result.tsv \
--contigs $FILEPATH/assembly/${j}/contigs.fasta \
--write_bins \
--write_bin_evals \
--threads ${SLURM_CPUS_PER_TASK} \
--outputbasename $FILEPATH/binning/${j}/mbin

## -------------------------------------------------------------- ##
## CHECKM2 - binning QC
## -------------------------------------------------------------- ##

# make output directory
mkdir $FILEPATH/binning_stats

checkm2 predict \
--threads ${SLURM_CPUS_PER_TASK} \
--extension .fa \
--database_path /path_to_database/uniref100.KO.1.dmnd \
--input $FILEPATH/binning/${j}/mbin_DASTool_bins \
--output-directory $FILEPATH/binning_stats/${j}

## ------------------------------------------------------------- ##
## THE END
## ------------------------------------------------------------- ##
