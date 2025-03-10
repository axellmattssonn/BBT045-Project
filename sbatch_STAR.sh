#!/usr/bin/env bash

## This specifies the project and what cluster to run on. Don't change it for this course!
#SBATCH -A C3SE2024-2-16 -p vera
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -t 13:00:00
#SBATCH -J STAR
#SBATCH --mail-user=axelmat@student.chalmers.se
#SBATCH --mail-type=ALL
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out

###
#
# Title: STAR.sh
# Date: 2024.02.23
# Author: Gastric Warriors inspired Johan Bengtsson-Palme and Vi Varga
#
# Description: 
# This script will run STAR alignment on a set of input samples in a loop.
#
# Usage: 
# sbatch sbatch_STAR.sh
#
###

### Set parameters
WORKDIR=/cephyr/NOBACKUP/groups/bbt045_2025/groups/group_gaswar/
INDEX=${WORKDIR}/human_index/
REF=${WORKDIR}/Reference/Homo_sapiens.GRCh38.113.gtf

# Container file
CONTAINER_LOC=${WORKDIR}/gaswar_container.sif

# Temporary working directory
WORKING_TMP=$TMPDIR/JOB_TMP

### Load modules
module purge
#module load MODULE_NAME/module.version ...;

### Prepare working directory on compute node
mkdir -p $WORKING_TMP
cd $WORKING_TMP

### Running STAR

## Identify a list of trimmed FASTQ forward read files
FILE_LIST=$(ls /cephyr/NOBACKUP/groups/bbt045_2025/groups/group_gaswar/trimmed_fastq/*_1_trimmed_fastq.gz)

for file1 in $FILE_LIST; do
    # Extract the basename (sample ID)
    file=$(basename "$file1" _1_trimmed_fastq.gz)
    
    # Define the paired-end file with full path
    file2="/cephyr/NOBACKUP/groups/bbt045_2025/groups/group_gaswar/trimmed_fastq/${file}_2_trimmed_fastq.gz"

    echo "Running STAR for: $file1 and $file2"

    # Run STAR inside the container
    apptainer exec $CONTAINER_LOC STAR \
        --runMode alignReads \
        --genomeDir $INDEX \
        --outSAMtype BAM SortedByCoordinate \
        --readFilesIn ${file1} ${file2} \
        --runThreadN 8 \
        --readFilesCommand gunzip -c \
        --outFileNamePrefix "${WORKING_TMP}/${file}_"

done

OUTPUT=${WORKDIR}/gene_counts.txt

apptainer exec $CONTAINER_LOC featureCounts -T 8 -p -t exon -g gene_id \
    -a $REF \
    -o $OUTPUT \
    $WORKING_TMP/*.bam

### Copy relevant files back, this is good practice but will actually not do anything for this specific script
cp -v $WORKING_TMP/*.bam $WORKDIR/Align/
cp -v $WORKING_TMP/*.out $WORKDIR/Align/
cp -v $WORKING_TMP/*.tab $WORKDIR/Align/