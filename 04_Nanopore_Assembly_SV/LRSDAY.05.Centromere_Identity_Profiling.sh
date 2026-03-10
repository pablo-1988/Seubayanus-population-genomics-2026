#!/bin/bash
# =============================================================================
# LRSDAY Step 05: Centromere Identity Profiling
# Long-Read Sequencing Data Analysis Workflow (LRSDAY) — S. eubayanus
# =============================================================================
# Description:
#   Identifies putative centromere sequences in the assembly by running
#   Exonerate to align known centromere query sequences against each
#   chromosome. The resulting alignments define centromere coordinates
#   for downstream annotation steps.
#
# Example strain: UCD646 (Holarctic-admixed S. eubayanus, Oxford Nanopore)
#
# Input:
#   Assembly FASTA + centromere query sequences (FASTA)
#
# Output:
#   Centromere coordinate table / alignment output file
#
# Dependencies:
#   Exonerate
#
# LRSDAY reference: https://github.com/yjx1217/LRSDAY
# =============================================================================
set -e -o pipefail
#######################################
# load environment variables for LRSDAY
source ./../../env.sh

#######################################
# set project-specific variables
prefix="UCD646" # The file name prefix for the processing sample. Default = "SK1" for the testing example.
genome="./../04.Reference-guided_Assembly_Scaffolding/$prefix.assembly.ref_based_scaffolded.fa" # The file path of the input genome assembly
query="./../00.Reference_Genome/CBS12357.fasta" # S. cerevisiae S288C reference centromere sequences based on Yue et al. (2017) Nature Genetics. 
debug="no" # Whether to keep intermediate files for debugging. Use "yes" if prefer to keep intermediate files, otherwise use "no". Default = "no".

######################################
# process the pipeline
$exonerate_dir/exonerate --showvulgar no --showcigar no --showalignment no --showtargetgff yes --bestn 1 $query $genome >$prefix.centromere.exonerate.gff
perl $LRSDAY_HOME/scripts/exonerate_gff2gff3.pl  -i $prefix.centromere.exonerate.gff -o $prefix.centromere.gff3.tmp -t $prefix
perl $LRSDAY_HOME/scripts/tidy_maker_gff3.pl -r $genome -i  $prefix.centromere.gff3.tmp -o  $prefix.centromere.gff3 -t $prefix

# clean up intermediate files
if [[ $debug == "no" ]]
then
    rm $prefix.centromere.exonerate.gff
    rm $prefix.centromere.gff3.tmp
fi

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
