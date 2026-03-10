#!/bin/bash
# =============================================================================
# LRSDAY Step 00: Prepare S. cerevisiae Reference Genome
# Long-Read Sequencing Data Analysis Workflow (LRSDAY) — S. eubayanus
# =============================================================================
# Description:
#   Copies and decompresses the S288C S. cerevisiae reference genome FASTA
#   for use in downstream scaffolding, annotation, and comparative genomics
#   steps throughout the LRSDAY workflow.
#
# Example strain: UCD646 (Holarctic-admixed S. eubayanus, Oxford Nanopore)
#
# Input:
#   Compressed reference genome FASTA (.fasta.gz or .fa.gz)
#
# Output:
#   Decompressed reference genome FASTA file
#
# Dependencies:
#   gzip
#
# LRSDAY reference: https://github.com/yjx1217/LRSDAY
# =============================================================================
set -e -o pipefail

#######################################
# load environment variables for LRSDAY
source ./../../env.sh

#######################################
# set project-specific variables

#######################################
# process the pipeline

cp $LRSDAY_HOME/data/S288C.ASM205763v1.fa.gz .
cp $LRSDAY_HOME/data/S288C.ASM205763v1.noncore_masked.fa.gz .
gunzip *.gz

############################
# checking bash exit status
if [[ $? -eq 0 ]]
then
    echo ""
    echo "LRSDAY message: This bash script has been successfully processed! :)"
    echo ""
    echo ""
    exit 0
fi
############################
