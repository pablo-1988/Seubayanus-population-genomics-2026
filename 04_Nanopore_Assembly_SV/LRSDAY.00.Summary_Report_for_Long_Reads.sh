#!/bin/bash
# =============================================================================
# LRSDAY Step 00: Summary Report for Long Reads
# Long-Read Sequencing Data Analysis Workflow (LRSDAY) — S. eubayanus
# =============================================================================
# Description:
#   Generates quality control statistics and plots for Nanopore long reads
#   using NanoPlot, including read length distribution, N50, and per-read
#   quality scores. The HTML report aids in assessing sequencing run quality.
#
# Example strain: UCD646 (Holarctic-admixed S. eubayanus, Oxford Nanopore)
#
# Input:
#   Nanopore FASTQ file(s) (raw or filtered)
#
# Output:
#   NanoPlot HTML report with QC statistics and visualizations
#
# Dependencies:
#   NanoPlot
#
# LRSDAY reference: https://github.com/yjx1217/LRSDAY
# =============================================================================
set -e -o pipefail
#######################################
# load environment variables for LRSDAY
source ./../../env.sh

#######################################
# set project-specific variables
long_reads_in_fastq="UCD646.porechop.fastq.gz" # The fastq file of long-reads (in fastq or fastq.gz format). Default = "SK1.filtered_subreads.fastq.gz"
prefix="UCD646" # The file name prefix for output files of the testing example. 
threads=1 # The number of threads to use. Default = "1".

#######################################
# process the pipeline

#source $nanoplot_dir/activate
$nanoplot_dir/NanoPlot \
    --threads $threads \
    --fastq $long_reads_in_fastq \
    --minlength 0 \
    --drop_outliers \
    --N50 \
    -o "${prefix}_Long_Reads_Summary_Report_out"

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
