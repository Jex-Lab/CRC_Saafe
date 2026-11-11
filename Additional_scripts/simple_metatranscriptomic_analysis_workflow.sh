#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=
#SBATCH --mem=
#SBATCH --time=
#SBATCH --array=0-99%10
#SBATCH --job-name=
#SBATCH --output=humann_%A-%a.out

# set path
FILEPATH=`pwd -P`

# read list of sample names
readarray -t samples < $FILEPATH/sample_list.txt

# assign each sample to an array task
j=${samples[${SLURM_ARRAY_TASK_ID}]}

## -------------------------------------------------------------- ##
## KNEADDATA - trim and remove contamination (human and rRNA)
## -------------------------------------------------------------- ##

# make output directory
mkdir -p $FILEPATH/clean/${j}

kneaddata \
--threads ${SLURM_CPUS_PER_TASK} \
--input1 $FILEPATH/data_raw/${j}_combined_R1.fastq.gz \
--input2 $FILEPATH/data_raw/${j}_combined_R2.fastq.gz \
--reference-db /data_reference/human_index/dna/bowtie2 \
--reference-db /data_reference/human_index/rna/bowtie2 \
--reference-db /data_reference/kneaddata_DBs/silva \
--output-prefix ${j}_clean \
--remove-intermediate-output \
--output $FILEPATH/clean/${j} \
--trimmomatic-options="ILLUMINACLIP:/data_reference/adapters/illuminaClipping.fa:2:40:20" \
--trimmomatic-options="LEADING:3" \
--trimmomatic-options="TRAILING:3" \
--trimmomatic-options="HEADCROP:14" \
--trimmomatic-options="SLIDINGWINDOW:4:15" \
--trimmomatic-options="MINLEN:25" \
--bowtie2-options="--very-sensitive"

# delete redundant working files
rm $FILEPATH/clean/${j}/${j}*contam*
rm $FILEPATH/clean/${j}/${j}*unmatc*

# compress output
pigz -p ${SLURM_CPUS_PER_TASK} $FILEPATH/clean/${j}/${j}_clean_paired_1.fastq
pigz -p ${SLURM_CPUS_PER_TASK} $FILEPATH/clean/${j}/${j}_clean_paired_2.fastq

## ------------------------------------------------------------- ##
## HUMANN4 - identify metabolic potential and activity
## ------------------------------------------------------------- ##

# make output directory
mkdir -p $FILEPATH/humann4_out/${j}

# concatenate fwd and rev reads
cat $FILEPATH/clean/${j}/${j}_clean_paired_1.fastq.gz $FILEPATH/clean/${j}/${j}_clean_paired_2.fastq.gz > $FILEPATH/humann4_out/${j}/${j}_clean_merged.fastq.gz

# calculate average read length
arl=`seqkit stats $FILEPATH/humann4_out/${j}/${j}_clean_merged.fastq.gz | awk '{print $7}' | tail -n 1`

humann4 \
--threads ${SLURM_CPUS_PER_TASK} \
--input $FILEPATH/humann4_out/${j}/${j}_clean_merged.fastq.gz \
--output $FILEPATH/humann4_out/${j} \
--bypass-nucleotide-search \
--average-read-length ${arl} \
--count-normalization Counts \
--resume

## ------------------------------------------------------------- ##
## THE END
## ------------------------------------------------------------- ##
