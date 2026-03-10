#!/bin/bash
# =============================================================================
# LRSDAY Step 02: Long-Read-Based Assembly Polishing
# Long-Read Sequencing Data Analysis Workflow (LRSDAY) — S. eubayanus
# =============================================================================
# Description:
#   Polishes the raw long-read assembly using the original Nanopore reads
#   through iterative rounds of Racon consensus correction followed by
#   Medaka neural-network-based polishing to reduce insertion/deletion errors.
#
# Example strain: UCD646 (Holarctic-admixed S. eubayanus, Oxford Nanopore)
#
# Input:
#   Raw genome assembly FASTA + filtered Nanopore FASTQ reads
#
# Output:
#   Long-read polished genome assembly FASTA
#
# Dependencies:
#   Racon, Medaka
#
# LRSDAY reference: https://github.com/yjx1217/LRSDAY
# =============================================================================
set -e -o pipefail
##########################################
# load environment variables for LRSDAY
source ./../../env.sh

###########################################
# set project-specific variables

prefix="UCD646" # The file name prefix for the processing sample. Please avoid the character '.' in prefix. Default = "SK1" for the testing example.
input_assembly="./../01.Long-read-based_Genome_Assembly/$prefix.assembly.raw.fa" # The file path of the input raw long-read-based assembly for polishing.
long_reads_in_fastq="./../00.Long_Reads/${prefix}.filtlong.fastq.gz" # The file path of the long-read fastq file. 

polisher="racon-medaka" # The long-read-based polisher to use: "quiver" (for PacBio RSII reads), "arrow" (for PacBio Sequel reads), "nanopolish" (for raw nanopore fast5 reads), "racon-medaka" (for basecalled nanopore fastq reads), or "marginpolish" (for basecalled nanopore fastq reads). Default = "quiver" for the testing example.
guppy_basecalling_model="r941_flip235" # The guppy basecalling model to use for medaka. Supported values include: "r941_min_fast" for guppy (version_number >= 3.0.3) in fast mode, "r941_min_high" for guppy (version_number >= 3.0.3) in high accuracy mode, "r941_flip235" for guppy (2.3.5 <= version_number <3.0.3), "r941_flip213" for guppy (2.1.3 <= version_number <2.3.5), and "r941_trans" for albacore or guppy (version_number < 2.1.3). This option is only needed when polisher="racon-medaka".
threads=4  # The number of threads to use. Default = "1".
ploidy=2 # The ploidy status of the sequenced genome. Use "1" for haploid genome and "2" for diploid genome. Currently not supported when "polisher="racon-medaka". Default = "1" for the testing example.
rounds_of_successive_polishing=3 # The number of total rounds of long-read-based assembly polishing. Default = "1" for the testing example.
debug="no" # Use "yes" if prefer to keep intermediate files, otherwise use "no". Default = "no"

###########################################
# process the pipeline

cp $input_assembly $prefix.assembly.tmp.fa

mkdir tmp

if [[ $polisher == "racon-medaka" ]]
then
    source $miniconda2_dir/activate $conda_medaka_dir/../../conda_medaka_env
    for i in $(seq 1 1 $rounds_of_successive_polishing)
    do
	$minimap2_dir/minimap2 -t $threads -ax map-ont $prefix.assembly.tmp.fa $long_reads_in_fastq > $prefix.minimap2.round_${i}.sam
	$racon_dir/racon -t $threads $long_reads_in_fastq $prefix.minimap2.round_${i}.sam $prefix.assembly.tmp.fa > $prefix.assembly.racon.round_${i}.fa
	if [[ $debug == "no" ]]
	then
	    rm $prefix.minimap2.round_${i}.sam
	fi
	rm $prefix.assembly.tmp.fa
	cp $prefix.assembly.racon.round_${i}.fa $prefix.assembly.tmp.fa
    done
    for i in $(seq 1 1 $rounds_of_successive_polishing)
    do
	$conda_medaka_dir/medaka_consensus -i $long_reads_in_fastq -d $prefix.assembly.tmp.fa -o ${prefix}_medaka_out_round_${i} -t $threads -m $guppy_basecalling_model
	rm $prefix.assembly.tmp.fa
	rm $prefix.assembly.tmp.fa.mmi
	rm $prefix.assembly.tmp.fa.fai
        perl $LRSDAY_HOME/scripts/tidy_fasta_for_medaka.pl -i ${prefix}_medaka_out_round_${i}/consensus.fasta -o $prefix.assembly.medaka.round_${i}.fa 
	cp $prefix.assembly.medaka.round_${i}.fa $prefix.assembly.tmp.fa
	if [[ $debug == "no" ]]
	then
	    rm -r ${prefix}_medaka_out_round_${i}
	fi
    done
    ln -s $prefix.assembly.medaka.round_${rounds_of_successive_polishing}.fa $prefix.assembly.long_read_polished.fa
    rm $prefix.assembly.tmp.fa
    source $miniconda2_dir/deactivate
fi

rm -r tmp

# clean up intermediate files
if [[ $debug == "no" ]]
then
    echo "clean up"
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
