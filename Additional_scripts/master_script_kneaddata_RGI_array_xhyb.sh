#!/bin/bash
#SBATCH --mem=150G
#SBATCH --cpus-per-task=40
#SBATCH --time=48:00:00
#SBATCH --array=1-4
#SBATCH --job-name=knead_rgi_sub
#SBATCH --output=/vast/scratch/users/azeez.i/20260513_rgi_test_backup6/rgi_test/card_DBs/test2/slurm_log/%x-%A-%a.out

# Lines starting with #SBATCH are for SLURM job management systems
# and may be removed if user is not submitting jobs to SLURM

####################################################################
#
# master_script_filenames_kneaddata_metaphlan.sh
# Created May 2026 by Imadh Abdul Azeez
#
####################################################################
#
# This script sets up and runs raw fastq sequences through the bioBakery kneadata and RGI pipelines for AMR analysis, specifically for job arrays
# Prior to this script, the raw fastq files have to be sorted by filenames into individual directories, using create_array_dirs_filenames.sh
#
# The steps are:
#   1.  Running kneaddata to obtain cleaned paired end reads
#   2.  Running RGI to obtain AMR profile
#
# NOTE: BEFORE running this script, please run create_array_dirs_filenames.sh
#
#######################################################################
#
echo -e "NOTE: Before running this script, please ensure filename array directoies have been made and that kneaddata and metaphlan has been installed with databases.\n\n"
#
# VARIABLES - set starting location and starting files location pathways
# Each Slurm array task chooses an input subdirectory based on its task ID
# NOTE: ensure that before creating array directories, that fastq files are gunzipped
INPUT_DIR=/vast/scratch/users/azeez.i/20260513_rgi_test_backup6/rgi_test/card_DBs/test2/input_file_filename                                           # MODIFY THIS IF NEEDED
ARRAY_SAMPLE_PREFIX=`ls $INPUT_DIR | sed ${SLURM_ARRAY_TASK_ID}'q;d'`
INPUT_DIR=$INPUT_DIR/$ARRAY_SAMPLE_PREFIX
OUT_DIR=/vast/scratch/users/azeez.i/20260513_rgi_test_backup6/rgi_test/card_DBs/test2/separated_outputs/output_$ARRAY_SAMPLE_PREFIX
mkdir -p $OUT_DIR

# Databases
knead_db=/stornext/General/data/academic/lab_jex/Imadh/20251218_from_vast_project_from_lab_share/databases/biobakery/kneaddata/kneaddata_human_genome_database
knead_db2=/stornext/General/data/academic/lab_jex/Imadh/20251218_from_vast_project_from_lab_share/databases/biobakery/kneaddata/kneaddata_human_transcriptome_database

# number of threads
threads=${SLURM_CPUS_PER_TASK} #`getconf _NPROCESSORS_ONLN`

STEP_1="$OUT_DIR/step_1_output"
STEP_2="$OUT_DIR/step_2_output"
STEP_3="$OUT_DIR/step_3_output"


mkdir -p ${STEP_1} ${STEP_2}

####################################################################
#STEP 0.1: create/read checkpoint

printf "\nStep 0.1: Checking for the presence of the checkpoint file.\n"
if [ ! -f "$INPUT_DIR/checkpoints" ]
  then
    printf "\tThe file 'checkpoints' does not exist in the input directory, creating...\n"
    touch "$INPUT_DIR/checkpoints"
else
    printf "\tThe file 'checkpoints' already exists in the input directory.\n"
fi

####################################################################
#
# STEP 1: RUNNING FASTQ FILES THROUGH KNEADDATA

Step=$(grep "KNEAD" $INPUT_DIR/checkpoints)
if [ "${Step}" != "KNEAD" ]
  then

$MKDIR $STEP_1
paired=false
for f in $INPUT_DIR/*R1*q
do
    f2=`echo $f | awk -F "R1" '{print $1 "R2" $2}'`
    if [ -f $f2 ]; then
      paired=true
      kneaddata --threads $threads \
        --input $f \
        --input $f2 \
        -db $knead_db \
        -db $knead_db2 \
        --output $INPUT_DIR \
        --trimmomatic-options="ILLUMINACLIP:/vast/projects/SAMSA/samsa2/programs/Trimmomatic-0.36/adapters/illuminaClipping.fa:2:40:10" \
        --bypass-trf \
        --trimmomatic-options="LEADING:3" \
        --trimmomatic-options="TRAILING:3" \
        --trimmomatic-options="SLIDINGWINDOW:4:15" \
        --trimmomatic-options="MINLEN:25"
    else
      kneaddata --threads $threads \
        --input $f \
        -db $knead_db \
        -db $knead_db2 \
        --output $INPUT_DIR \
        --trimmomatic-options="ILLUMINACLIP:/vast/projects/SAMSA/samsa2/programs/Trimmomatic-0.36/adapters/illuminaClipping.fa:2:40:10" \
        --bypass-trf \
        --trimmomatic-options="LEADING:3" \
        --trimmomatic-options="TRAILING:3" \
        --trimmomatic-options="SLIDINGWINDOW:4:15" \
        --trimmomatic-options="MINLEN:25"
    fi
done

if $paired; then
  mv $INPUT_DIR/*"_paired_1.fastq"* $STEP_1
  mv $INPUT_DIR/*"_paired_2.fastq"* $STEP_1
else
  mv $INPUT_DIR/*"kneaddata.fastq" $STEP_1
fi


printf "KNEAD\n" >>$INPUT_DIR/checkpoints

else
  printf  "\tThe variable KNEAD is in the checkpoint file. STEP 1 will be skipped.\n"
fi


####################################################################
#
# STEP 2: RUNNING CLEAG FASTQ FILES THROUGH RGI
# 

# Set temporary directory in VAST scratch space to deal with larger samples
mkdir -p /vast/scratch/users/azeez.i/tmp/$ARRAY_SAMPLE_PREFIX
export TMPDIR=/vast/scratch/users/azeez.i/tmp/$ARRAY_SAMPLE_PREFIX
echo $TMPDIR

Step=$(grep "RGI" $INPUT_DIR/checkpoints)
if [ "${Step}" != "RGI" ]
  then

  $MKDIR $STEP_2
  if $paired; then
    for file in $STEP_1/*_paired_1.fastq
    do
      f2=`echo $file | awk -F "paired_1.fastq" '{print $1 "paired_2.fastq"}'`
      rgi bwt -1 $file -2 $f2 --output $STEP_2/$ARRAY_SAMPLE_PREFIX -n $threads --mapq 30 --coverage 0.8 --local -n 16
    done
  else
    for file in $INPUT_DIR/*kneaddata.fastq
    do
      rgi bwt -1 $file --output $STEP_2/$ARRAY_SAMPLE_PREFIX --local
    done
  fi

  printf "RGI\n" >>$INPUT_DIR/checkpoints

else
  printf  "\tThe variable RGI is in the checkpoint file. STEP 2 will be skipped.\n"
fi


####################################################################
echo "Master bash script finished running at: "; date
exit 0
####################################################################
