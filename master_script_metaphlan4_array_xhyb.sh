#!/bin/bash
#SBATCH --mem=150G
#SBATCH --cpus-per-task=40
#SBATCH --time=48:00:00
#SBATCH --array=1-32%12
#SBATCH --job-name=xhyb_mpa
#SBATCH --output=/vast/scratch/users/azeez.i/20260513_rgi_test_backup6/rgi_test/card_DBs/test_for_prep_in_run002/slurm_logs_subsample_spike/%x-%A-%a.out

# Lines starting with #SBATCH are for SLURM job management systems
# and may be removed if user is not submitting jobs to SLURM

####################################################################
#
# master_script_metaphlan4_array_xhyb.sh
# Created May 2026 by Imadh Abdul Azeez
#
####################################################################
#
# This is the second script to after the first script (master_script_kneaddata_RGI_array_xhyb.sh) sets up and runs raw fastq sequences through the bioBakery kneadata and RGI pipelines
# This script is applied afterwards for for microbiome community analysis, specifically for job arrays
# Prior to this script, the raw fastq files have to be sorted by filenames into individual directories, using create_array_dirs_filenames.sh, then the master_script_kneaddata_RGI_array_xhyb.sh and finally this script
#
# The steps in this script are:
#   1.  Running Metaphlan to obtain community profile
#
# NOTE: BEFORE running this script, please run create_array_dirs_filenames.sh and master_script_kneaddata_RGI_array_xhyb.sh
# NOTE: this script will output the community profile in a step_3_output sub-directory in each output task ID directory in separated outputs
#
#######################################################################
#
echo -e "NOTE: Before running this script, please ensure filename array directoies have been made and that kneaddata and metaphlan has been installed with databases.\n\n"
#
# VARIABLES - set starting location and starting files location pathways
# Each Slurm array task chooses an input subdirectory based on its task ID
# NOTE: ensure that before creating array directories, that fastq files are gunzipped
INPUT_DIR=/vast/scratch/users/azeez.i/20260513_rgi_test_backup6/rgi_test/card_DBs/test_for_prep_in_run002/input_files_filenames_subsample_spike                                           # MODIFY THIS IF NEEDED
ARRAY_SAMPLE_PREFIX=`ls $INPUT_DIR | sed ${SLURM_ARRAY_TASK_ID}'q;d'`
INPUT_DIR=$INPUT_DIR/$ARRAY_SAMPLE_PREFIX
OUT_DIR=/vast/scratch/users/azeez.i/20260513_rgi_test_backup6/rgi_test/card_DBs/test_for_prep_in_run002/seperated_outputs_test_subsample_spike/output_$ARRAY_SAMPLE_PREFIX
# mkdir -p $OUT_DIR

# Databases
# knead_db=/stornext/General/data/academic/lab_jex/Imadh/20251218_from_vast_project_from_lab_share/databases/biobakery/kneaddata/kneaddata_human_genome_database
# knead_db2=/stornext/General/data/academic/lab_jex/Imadh/20251218_from_vast_project_from_lab_share/databases/biobakery/kneaddata/kneaddata_human_transcriptome_database

# number of threads
threads=${SLURM_CPUS_PER_TASK} #`getconf _NPROCESSORS_ONLN`

STEP_1="$OUT_DIR/step_1_output"
STEP_2="$OUT_DIR/step_2_output"
STEP_3="$OUT_DIR/step_3_output"


mkdir -p ${STEP_3}

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
# STEP 3: RUNNING CLEAG FASTQ FILES THROUGH METAPHLAN
#

# Since this script is expected to use paired input files, set paired to true. Modify script if we want to modularly choose between paired and single-end
paired=true

Step=$(grep "META" $INPUT_DIR/checkpoints)
if [ "${Step}" != "META" ]
  then

$MKDIR $STEP_3
if $paired; then
  for file in $STEP_1/*_paired_1.fastq
  do
    f2=`echo $file | awk -F "paired_1.fastq" '{print $1 "paired_2.fastq"}'`
    metaphlan $file,$f2 --input_type fastq --bowtie2db /vast/scratch/users/azeez.i/20260519_humann4_update_db_for_strat_backup3/humann4_update_db_for_strat/metaphlan_databases/vOct22 --index mpa_vOct22_CHOCOPhlAnSGB_202403 --bowtie2out $STEP_3/$ARRAY_SAMPLE_PREFIX"_metagenome.bz2" -o $STEP_3/$ARRAY_SAMPLE_PREFIX"_profiled_metagenome.txt"
  done
else
  for file in $STEP_1/*kneaddata.fastq
  do
    new_name=`echo $file | awk -F "cleaned" '{print $1 "profiled_metagenome.txt"}'`
    metaphlan $file --input_type fastq  $STEP_3/${new_name##*/}
  done
fi

printf "META\n" >>$INPUT_DIR/checkpoints

else
  printf  "\tThe variable META is in the checkpoint file. STEP 2 will be skipped.\n"
fi

####################################################################
echo "Master bash script finished running at: "; date
exit 0
####################################################################
