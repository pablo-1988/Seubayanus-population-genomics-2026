#!/bin/bash
# =============================================================================
# LRSDAY Step 07 (Part 1): Supervised Final Assembly — Generate Modification List
# Long-Read Sequencing Data Analysis Workflow (LRSDAY) — S. eubayanus
# =============================================================================
# Description:
#   Generates a modification list template by aligning the scaffolded assembly
#   against the reference genome with MUMmer. The template guides manual
#   curation by listing contigs/scaffolds that may need relabeling, reordering,
#   or splitting before producing the final assembly.
#
# Example strain: UCD646 (Holarctic-admixed S. eubayanus, Oxford Nanopore)
#
# Input:
#   Scaffolded assembly FASTA
#
# Output:
#   Modification list template file for manual review and editing
#
# Dependencies:
#   MUMmer, custom LRSDAY scripts
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
genome="./../06.Mitochondrial_Genome_Assembly_Improvement/$prefix.assembly.mt_improved.fa" # The file name of the input genome assembly.


#######################################
# process the pipeline
# Step 1:
echo "#original_name,orientation,new_name" > ${prefix}.assembly.modification.list
cat $genome |egrep ">"|sed "s/>//gi"|awk '{print $1 ",+," $1}' >>${prefix}.assembly.modification.list

echo "################################"
echo "running LRSDAY.06.Supervised_Final_Assembly.1.sh > Done!"
echo "Please manually edit the generated $prefix.modification.list for relabeling/reordering contigs when necessary"
echo "Once you finish the editing, plase run the script LRSDAY.06.Supervised_Final_Assembly.2.sh."
echo "################################"

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
