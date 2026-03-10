#!/bin/bash
# =============================================================================
# LRSDAY Step 08: Centromere Annotation
# Long-Read Sequencing Data Analysis Workflow (LRSDAY) — S. eubayanus
# =============================================================================
# Description:
#   Re-annotates centromere sequences in the final curated assembly by running
#   Exonerate with centromere query sequences. The resulting coordinates are
#   formatted as BED and GFF3 files for integration into the final annotation.
#
# Example strain: UCD646 (Holarctic-admixed S. eubayanus, Oxford Nanopore)
#
# Input:
#   Final assembly FASTA + centromere query sequences (FASTA)
#
# Output:
#   Centromere annotation in BED and/or GFF3 format
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
genome="./../07.Supervised_Final_Assembly/$prefix.assembly.final.fa" # The path of the input genome assembly.
query="./../00.Reference_Genome/CBS12357.fasta" # The S. cerevisiae S288C reference centromere sequences based on Yue et al. (2017) Nature Genetics.
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
