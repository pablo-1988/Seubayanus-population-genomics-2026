#!/bin/bash
# =============================================================================
# TreeMix Pipeline — Population Graph Analysis of S. eubayanus
# =============================================================================
# Description:
#   Runs a complete TreeMix analysis starting from a quality-filtered VCF file.
#   Steps include:
#     1. Subsetting pure-lineage strains from the VCF
#     2. Filtering for biallelic sites only
#     3. Renaming chromosomes (required for LD pruning with PLINK)
#     4. Linkage disequilibrium (LD) pruning
#     5. Converting the pruned VCF to TreeMix input format
#     6. Running TreeMix with 0-8 migration events (10 iterations each)
#     7. Generating a consensus tree using PHYLIP Consense
#
# Input:
#   - hqvariants_yayaout.recode.vcf : quality-filtered VCF (no indels)
#   - keep.txt        : list of sample names to retain (pure strains)
#   - treeorder.clust : 3-column tab-separated file (sample, sample, population)
#                       used to define population groups for TreeMix
#
# Output:
#   - treemix.{m}.{iter} files : TreeMix results for m migrations, iteration iter
#   - Consensus tree files per migration value (via PHYLIP Consense)
#
# Dependencies:
#   vcftools, PLINK, TreeMix (>=1.13), PHYLIP (consense), bgzip, tabix, bcftools
#   Conda environments: base (vcftools/plink), treemix (separate env)
#
# Usage:
#   bash pipeline_treemix.sh
#   Note: Update file paths and chromosome names (sed commands) as needed
# =============================================================================

# --- Step 1: Subset samples to pure strains only ---
# keep.txt contains the list of sample names to retain
vcftools --vcf hqvariants_yayaout.recode.vcf --keep keep.txt --recode --recode-INFO-all --non-ref-ac-any 1 --out subset.vcf

# --- Step 2: Retain only biallelic SNP sites (no multiallelic, no missing data) ---
vcftools --max-missing 1 --max-alleles 2 --vcf subset.vcf.recode.vcf --recode --recode-INFO-all --out biallelicvariants

# --- Step 3: Rename chromosome identifiers to numbers ---
# PLINK requires numeric chromosome names for LD pruning
# Repeat this sed command for each chromosome, changing the pattern accordingly
sed -i 's/YALI0A/1/g' biallelicvariants.recode.vcf

# --- Step 4: LD pruning ---
# Download and run the LD pruning script (removes correlated SNPs)
wget https://github.com/joanam/scripts/raw/master/ldPruning.sh
chmod +x ldPruning.sh
./ldPruning.sh biallelicvariants.recode.vcf

# Decompress the LD-pruned VCF
gzip -d biallelicvariants.recode.LDpruned.vcf.gz

# --- Step 5: Convert VCF to TreeMix input format ---
# treeorder.clust: tab-separated file with 3 columns (sample, sample, population)
./vcf2treemix.sh biallelicvariants.recode.LDpruned.vcf treeorder.clust

# --- Step 6: Run TreeMix with 0–8 migration events, 10 iterations each ---
# Switch to treemix conda environment (incompatible with vcftools/plink env)
conda deactivate
conda activate treemix

# Loop over migration values (m) and iterations
# A different random seed is used per iteration (required for OptM analysis in R)
# Output files: treemix.{m}.{iter}
for i in {0..8}; do
  for iter in (1..10); do
    SEED=$RANDOM
    treemix -i biallelicvariants.recode.LDpruned.treemix.frq.gz \
            -m ${i} \
            -o treemix.${i}.${iter} \
            -bootstrap -k 500 -noss -seed $SEED \
            > treemix_${i}.${iter}_log &
  done
done
wait  # Wait for all background TreeMix processes to finish

# Note: After selecting the optimal m value using OptM (R package),
# rerun TreeMix 100 times at the optimal m to obtain the final topology.

# --- Step 7: Generate consensus trees using PHYLIP Consense ---
# Install PHYLIP if not already available
conda install phylip

# Organize TreeMix output by migration value (example for m=1)
mkdir m1   # Create directory for m=1 results
mv treemix.1.* m1/         # Move all m=1 output files
gzip -d m1/*.treeout.gz    # Decompress Newick tree files

# Concatenate tree lines into a single file for Consense input
# Consense requires an 'intree' file with all Newick trees concatenated
for file in m1/*.treeout; do head -n 1 "$file"; done > m1/atree.intree

# Repeat for other m values as needed (e.g., m3 folder)

# =============================================================================
# ADDITIONAL VCF PREPARATION STEPS
# (Used to annotate and reformat the VCF for TreeMix compatibility)
# =============================================================================

# Check allele frequencies in the VCF
vcftools --vcf treemix_vcf_f1.recode.vcf --out check --freq

# Decompress and recompress annotation table for tabix indexing
bgzip -d annotation_file.tab.gz
bgzip -c annotation_file.tab > annotation_file.tab.gz
tabix -s 1 -b 2 -e 2 annotation_file.tab.gz

# Annotate VCF with chromosome/position information
bcftools annotate -a annotation_file.tab.gz -c CHROM,POS treemix_vcf_f1.recode.vcf -o annotated_vcf.vcf

# Rename chromosome labels in VCF to numeric format required by PLINK/TreeMix
# Note: update the path and chromosome pattern as needed for your dataset
sed -i '' 's|CBS12357_Chr01_polished|1|g' /Users/pablovillarrealdiaz/Dropbox/the_vcf/treemix/biallelicvariants.recode.vcf
