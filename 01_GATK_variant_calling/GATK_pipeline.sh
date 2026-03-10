#!/bin/bash
# =============================================================================
# Population Genomics Pipeline — S. eubayanus
# GATK4 Variant Calling, Filtering, and Downstream Population Analyses
# =============================================================================
# Description:
#   Complete pipeline for SNP discovery and population genomics of S. eubayanus
#   strains sequenced by Illumina short-read technology. The pipeline follows
#   the GATK4 Best Practices workflow and includes:
#
#   PART 1 — Variant Calling:
#     - AddOrReplaceReadGroups (Picard) on aligned BAM files
#     - Per-sample variant calling (HaplotypeCaller in GVCF mode)
#     - Multi-sample genotyping (GenomicsDBImport + GenotypeGVCFs)
#     - VCF merging and quality filtering (vcftools hard filters)
#
#   PART 2 — Phylogenetic Analysis:
#     - VCF → PHYLIP format conversion (vcf2phylip)
#     - Maximum likelihood tree construction (IQTree, GTR+ASC model)
#
#   PART 3 — Population Structure:
#     - LD pruning with PLINK
#     - STRUCTURE (Bayesian clustering, K=2–10)
#     - ADMIXTURE (Maximum likelihood clustering, K=2–20)
#     - fineSTRUCTURE (chromosome painting-based fine-scale structure)
#     - fastStructure (variational Bayes, K=2–10, via Docker)
#
#   PART 4 — Gene Flow Analysis:
#     - TreeMix (population graph with migration edges, via STACKS/vcf conversion)
#
#   PART 5 — Gene Orthology:
#     - OrthoFinder (all-vs-all ortholog identification from GFF3 assemblies)
#
# Input:
#   - Sorted BAM files for each strain (from BWA-MEM alignment)
#   - Reference genome FASTA (S. eubayanus CBS12357 polished assembly)
#   - strains.txt : list of strain names
#   - mainparams / extraparams : STRUCTURE parameter files
#   - popstacks.txt : population mapping file for STACKS/TreeMix
#
# Output:
#   - all_raw_variants.vcf : raw multi-sample VCF
#   - filtered_snps.vcf : quality-filtered SNP VCF
#   - IQTree phylogenetic tree files (.treefile, .contree)
#   - ADMIXTURE .Q and .P files per K value
#   - fineSTRUCTURE XML results
#   - TreeMix output files
#   - OrthoFinder results folder
#
# Dependencies:
#   GATK4 (>=4.1.9), Picard, vcftools, bgzip, tabix, bcftools,
#   IQTree2, PLINK (>=1.9), ADMIXTURE, STRUCTURE, fineSTRUCTURE4,
#   fastStructure (Docker), TreeMix, STACKS (populations), OrthoFinder,
#   gffread, Python (vcf2phylip.py)
#
# Usage:
#   Run each section sequentially. Update file paths and parameters
#   (chromosome names, sample counts, K values) as needed for your dataset.
# =============================================================================

# =============================================================================
# PART 1: VARIANT CALLING
# =============================================================================

# Build sequence dictionary for the reference FASTA (required by GATK)
# Replace with your reference FASTA path
gatk CreateSequenceDictionary -R lach.ci35460_ref.fa


# Add read group information to BAM files (required by GATK HaplotypeCaller)
# Repeat for each strain, replacing CL1105 with the actual strain name
java -jar build/libs/picard.jar AddOrReplaceReadGroups -I ../Sequence/CL1105/Aln_1105.sorted.bam -O ../Sequence/CL1105/Aln_1105_sortedW@RG.bam -SORT_ORDER coordinate -RGID DontKnow -RGLB 01 -RGPL illumina -RGSM Sample -RGPU 001 -CREATE_INDEX True



# Step 1: Per-sample variant calling using HaplotypeCaller in GVCF mode
# Run for each strain and each chromosome separately (parallelizable)
# Example: chromosome LACI0A
for SAMPLE in ${cat strains}
do ../../gatk-4.1.9.0/gatk HaplotypeCaller --intervals LACI0A --emit-ref-confidence GVCF -R ../ref_genoma/lach.ci35450_ref.fa  -I ${SAMPLE}.bam -O ${SAMPLE}.vcf
done

# Step 2: Build variant database per chromosome (multi-sample consolidation)
# GenomicsDBImport combines per-sample GVCFs into a per-chromosome database

java -jar ../../gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar GenomicsDBImport --genomicsdb-workspace-path all_variants_chrA/ --intervals LACI0A   --arguments_file listLACI0A &
java -jar ../../gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar GenomicsDBImport --genomicsdb-workspace-path all_variants_chrB/ --intervals LACI0B   --arguments_file listLACI0B &
java -jar ../../gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar GenomicsDBImport --genomicsdb-workspace-path all_variants_chrC/ --intervals LACI0C   --arguments_file listLACI0C &
java -jar ../../gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar GenomicsDBImport --genomicsdb-workspace-path all_variants_chrD/ --intervals LACI0D   --arguments_file listLACI0D &
java -jar ../../gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar GenomicsDBImport --genomicsdb-workspace-path all_variants_chrE/ --intervals LACI0E   --arguments_file listLACI0E &
java -jar ../../gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar GenomicsDBImport --genomicsdb-workspace-path all_variants_chrF/ --intervals LACI0F   --arguments_file listLACI0F &
java -jar ../../gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar GenomicsDBImport --genomicsdb-workspace-path all_variants_chrG/ --intervals LACI0G   --arguments_file listLACI0G &
java -jar ../../gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar GenomicsDBImport --genomicsdb-workspace-path all_variants_chrH/ --intervals LACI0H   --arguments_file listLACI0H &

###../../gatk-4.1.8.0/gatk GenomicsDBImport -V NS79_chr03.g.vcf  --genomicsdb-update-workspace-path all_variants_chrB&#########

../../gatk-4.1.8.0/gatk GenomicsDBImport -V NS79_LACI0A.vcf  --genomicsdb-update-workspace-path all_variants_chrA &
../../gatk-4.1.8.0/gatk GenomicsDBImport -V NS79_LACI0A.vcf  --genomicsdb-update-workspace-path all_variants_chrB&
../../gatk-4.1.8.0/gatk GenomicsDBImport -V NS79_LACI0A.vcf  --genomicsdb-update-workspace-path all_variants_chrC&
../../gatk-4.1.8.0/gatk GenomicsDBImport -V NS79_LACI0A.vcf  --genomicsdb-update-workspace-path all_variants_chrD&
../../gatk-4.1.8.0/gatk GenomicsDBImport -V NS79_LACI0A.vcf  --genomicsdb-update-workspace-path all_variants_chrE&
../../gatk-4.1.8.0/gatk GenomicsDBImport -V NS79_LACI0A.vcf  --genomicsdb-update-workspace-path all_variants_chrF&
../../gatk-4.1.8.0/gatk GenomicsDBImport -V NS79_LACI0A.vcf  --genomicsdb-update-workspace-path all_variants_chrG&
../../gatk-4.1.8.0/gatk GenomicsDBImport -V NS79_LACI0A.vcf  --genomicsdb-update-workspace-path all_variants_chrH&



# Step 3: Joint genotyping — call genotypes from the multi-sample database


 ../gatk-4.1.9.0/gatk GenotypeGVCFs -R ../ref_genoma/lach.ci35450_ref.fa -V gendb://all_variants_chrA/  -G StandardAnnotation -O variants_chrA.vcf 2>chrA_log&
 ../gatk-4.1.9.0/gatk GenotypeGVCFs -R ../ref_genoma/lach.ci35450_ref.fa -V gendb://all_variants_chrB/  -G StandardAnnotation -O variants_chrB.vcf 2>chrB_log&
 ../gatk-4.1.9.0/gatk GenotypeGVCFs -R ../ref_genoma/lach.ci35450_ref.fa -V gendb://all_variants_chrC/  -G StandardAnnotation -O variants_chrC.vcf 2>chrC_log&
 ../gatk-4.1.9.0/gatk GenotypeGVCFs -R ../ref_genoma/lach.ci35450_ref.fa -V gendb://all_variants_chrD/  -G StandardAnnotation -O variants_chrD.vcf 2>chrD_log&
 ../gatk-4.1.9.0/gatk GenotypeGVCFs -R ../ref_genoma/lach.ci35450_ref.fa -V gendb://all_variants_chrE/  -G StandardAnnotation -O variants_chrE.vcf 2>chrE_log&
 ../gatk-4.1.9.0/gatk GenotypeGVCFs -R ../ref_genoma/lach.ci35450_ref.fa -V gendb://all_variants_chrF/  -G StandardAnnotation -O variants_chrF.vcf 2>chrF_log&
 ../gatk-4.1.9.0/gatk GenotypeGVCFs -R ../ref_genoma/lach.ci35450_ref.fa -V gendb://all_variants_chrG/  -G StandardAnnotation -O variants_chrG.vcf 2>chrG_log&
 ../gatk-4.1.9.0/gatk GenotypeGVCFs -R ../ref_genoma/lach.ci35450_ref.fa -V gendb://all_variants_chrH/  -G StandardAnnotation -O variants_chrH.vcf 2>chrH_log&

# Step 4: Merge per-chromosome VCFs into a single multi-sample VCF 

java -jar ../../gatk-4.1.9.0/gatk MergeVcfs -I variants_chrA.vcf -I variants_chrB.vcf -I variants_chrC.vcf -I variants_chrD.vcf -I variants_chrE.vcf -I variants_chrF.vcf -I variants_chrG.vcf -I variants_chrH.vcf -O all_raw_variants.vcf

##>quality "hard filters" for SNPs and INDELs

##../gatk-4.1.9.0/gatk SelectVariants -R ../ref_genoma/lach.ci35450_ref.fa -V all_raw_variants.vcf --select-type-to-include SNP -O raw_snps.vcf

#../gatk-4.1.9.0/gatk VariantFiltration -R ../ref_genoma/lach.ci35450_ref.fa -V raw_snps.vcf --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "Se_snp_filter" -O filtered_snps.vcf

###vcf-concensus
bgzip -c file.vcf > file.vcf.gz
tabix -p vcf file.vcf.gz

cat ref.fa | vcf-consensus file.vcf.gz > out.fa

#vcftools mean qual

vcftools --vcf b1_raw_variants.vcf --recode --recode-INFO-all --minQ 30 --min-meanDP 10 --remove-indels --out b1_snp

#First vcftools is used to filter SNPs, individuals, etc. For the first analyisis we will use SNPs that were called in all individuals, and we only retain biallelic sites

#Create a text file with the name of the individual to remove (remove)
vcftools --remove remove --vcf FILE.vcf --recode --recode-INFO-all --non-ref-ac-any 1 --out some_cidri

 
##We will generate a robust phylogenetic tree using Maximum Likehood method implemented in the linux program IQTREE. 
##First we transform the VCF file to phylip format using a
 custom python script developed by edgardomortiz, vcf2phylip, and then we use it for IQTREE

conda install -c bioconda iqtree
git clone https://github.com/edgardomortiz/vcf2phylip.git
#First we transform the VCF file to phylip format using a custom python script developed by edgardomortiz, and then we use it for IQTREE
python vcf2phylip.py -i filtered_snps_max.recode.vcf
iqtree -s filtered_snps_max.recode.min4.phy -st DNA -o CL1105.1 -m GTR+ASC -nt 8 #outgroup is CBS707 as is a different species
iqtree -s OUTFILE.min4.recode.min4.phy.varsites.phy -st DNA -o CBS707 -m GTR+ASC -nt 8 -bb 1000 -redo 


# =============================================================================
# PART 3: POPULATION STRUCTURE ANALYSIS
# =============================================================================

# --- STRUCTURE / ADMIXTURE ---
# Bayesian and maximum-likelihood clustering to infer ancestral population proportions
# Requires LD-pruned subset of SNPs (~10k) for computational efficiency 
#Structure is recommended to work with a subset of SNPs (~10k snps), which we can obtain by filtering out SNPs that are in linkage disequilibrium (LD) using the software PLINK. 
#As for determining population structure we do not need an outgroup, we can remove it using vcftools (preferentially using the first unfiltered vcf).

>Create a text file with the name of the individual to remove (remove)
vcftools --remove remove --vcf FILE.vcf --recode --recode-INFO-all --non-ref-ac-any 1 --out onlycidri

#Now to use PLINK, we modify the VCF file to write names to each SNP (requirement if we want to use the list that PLINK produces with the SNPs that are in LD) 
#We will use the following script available in gist: https://gist.github.com/janxkoci/25d495e6cb9f21d5ee4af3005fb3c77a#file-plink_pruning_prep-sh


#conda install tabix bcftools plink

./plink_pruning_prep.sh onlycidri_fil.recode.vcf
plink --vcf onlycidri.recode_annot.vcf --double-id --allow-extra-chr --indep-pairwise 50 10 0.2 --out onlycidri_fil_ldfilter

#Check if the prune.in file contain SNP ids#
head onlycidri_fil_ldfilter.prune.in

#Now we use this file to retain these positions from the vcf with vcftools. 
#However we need to change the underscore character that is located prior to the SNP position to a tab character For example: CBS12357_Chr12_polished_8383542 --> CBS12357_Chr12_polished 8383542

sed 's/_/\t/' onlycidri_fil_ldfilter.prune.in > onlycidri_fil_ldfilter.prune.in.vcftools
vcftools --vcf onlycidri_fil.recode.vcf --positions onlycidri_fil_ldfilter.prune.in.vcftools --recode-INFO-all --recode --out onlycidri_fil_ldfilter

#If we still have more than 20k SNPs, we can filter SNPs by distance using vcftools --thin option. Here we will filter out SNPs that are closer than 500 bp

vcftools --vcf onlycidri_fil_ldfilter.recode.vcf --thin 500 --recode --recode-INFO-all --out only_fil_ldfilter_thinned

#with PGDSpider transfor vcf to STRUCTURE format




for m in {2..10}
   do
   structure -m mainparams -p extraparams -K $m -L 5033 -N 354 -i plink.recode.strct_in -o structure_eub_K.${m}.rep1 &> log_K.${m}.rep1&
   sleep 10
   structure -m mainparams -p extraparams -K $m -L 5033 -N 354 -i plink.recode.strct_in -o structure_eub_K.${m}.rep2 &> log_K.${m}.rep2&
   sleep 10
   structure -m mainparams -p extraparams -K $m -L 5033 -N 354 -i plink.recode.strct_in -o structure_eub_K.${m}.rep3 &> log_K.${m}.rep3&
   sleep 10
   structure -m mainparams -p extraparams -K $m -L 5033 -N 354 -i plink.recode.strct_in -o structure_eub_K.${m}.rep4 &> log_K.${m}.rep4&
   sleep 10
   structure -m mainparams -p extraparams -K $m -L 5033 -N 354 -i plink.recode.strct_in -o structure_eub_K.${m}.rep5 &> log_K.${m}.rep5
   sleep 10  
   done 


##ADMIXTURE

plink --vcf  --double-id --allow-extra-chr --indep-pairwise 50 5 0.2 --maf 0.05 --out onlyeubayanus_ld --make-bed --threads 2

awk '{$1=0;print $0}' onlyeubayanus_ld.bim > onlyeubayanus_ld.bim.tmp
mv  onlyeubayanus_ld.bim.tmp  onlyeubayanus_ld.bim

for k in {16..20}
  do
  admixture -j5 --cv oporto_ld.bed $k > oporto_ld.$k.log
done



for k in `seq 1 10`; do for i in `seq 1 10`; do cp onlyeubayanus_ld.bed onlyeubayanus_ld.run_$i.bed && admixture --cv -j5 -s $RANDOM onlyeubayanus_ld.run_$i.bed $k \
| tee onlyeubayanus_ld.run_$i.$k.log && rm onlyeubayanus_ld.run_$i.bed; done; done


#Run Admixture from 3 to 10 pop
for i in {2..10}
do
 admixture --cv $FILE.bed $i > log${i}.out
done


#To identify the best value of k clusters which is the value with lowest cross-validation error, we need to collect the cv errors.
awk '/CV/ {print $3,$4}' *out | cut -c 4,7-20 > $FILE.cv.error
#o 
grep "CV" *out | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//'  > $FILE.cv.error

#To make plotting easier, we can make a file with the individual names in one column and the species names in the second column. 
awk '{split($1,name,"."); print $1,name[2]}' $FILE.nosex > $FILE.list

#Now we are ready to plot the results in R. To make it a bit easier, Joana Meier has written an R script for you that generates the plot. 
#It requires four arguments, the prefix for the ADMIXTURE output files (-p ), the file with the species information (-i ), 
#the maximum number of K to be plotted (-k 5), and a list with the populations or species separated by commas (-l <pop1,pop2...>). 
#The list of populations provided with -l gives the order in which the populations or species shall be plotted. Note, that alternatively, 
#if working with your own data, you could also try [this](https://github.com/ramachandran-lab/pong/blob/master/pong-manual.pdf) for plotting.

#You can get it with wget:
wget https://github.com/speciationgenomics/scripts/raw/master/plotADMIXTURE.r
chmod +x plotADMIXTURE.r

#Now, let’s run it like so:


Rscript plotADMIXTURE.r -p $FILE -i $FILE.list -k 5 -l #list




# --- fineSTRUCTURE ---
# Fine-scale population structure using chromosome painting (ChromoPainter + fineSTRUCTURE)

conda install snpsift plink
conda install -c compbiocore perl-switch
wget https://people.maths.bris.ac.uk/~madjl/finestructure/fs_4.1.1.zip
wget https://github.com/gusevlab/germline/raw/master/phasing_pipeline.tar.gz <- MAKE ALL inside folder
wget https://faculty.washington.edu/browning/beagle/recent.versions/beagle_3.0.4_05May09.zip <- JUST NEED TO COPY BEAGLE.JAR TO THE PHASING PIPELINE FOLDER



#Split VCF file by chromosome using Sift

SnpSift split onlycidri.recode.vcf

#Convert VCF files to PLINK .ped / .map format

plink --vcf onlycidri.recode.LACI0A.vcf --recode12 --allow-extra-chr --double-id --geno 1 --out LACI0A
plink --vcf onlycidri.recode.LACI0B.vcf --recode12 --allow-extra-chr --double-id --geno 1 --out LACI0B
plink --vcf onlycidri.recode.LACI0C.vcf --recode12 --allow-extra-chr --double-id --geno 1 --out LACI0C
plink --vcf onlycidri.recode.LACI0D.vcf --recode12 --allow-extra-chr --double-id --geno 1 --out LACI0D
plink --vcf onlycidri.recode.LACI0E.vcf --recode12 --allow-extra-chr --double-id --geno 1 --out LACI0E
plink --vcf onlycidri.recode.LACI0F.vcf --recode12 --allow-extra-chr --double-id --geno 1 --out LACI0F
plink --vcf onlycidri.recode.LACI0G.vcf --recode12 --allow-extra-chr --double-id --geno 1 --out LACI0G
plink --vcf onlycidri.recode.LACI0H.vcf --recode12 --allow-extra-chr --double-id --geno 1 --out LACI0H

#Convert files to chromopainter format
# fuera del ambuente conda (conda deactivate)

for i in {A B C D E F G H};do perl ./fs_4.1.1/plink2chromopainter.pl -p=LACI0${i}.ped -m=LACI0${i}.map -o=LACI0${i}.chromopainter -f;done

# entrar nuevamente en conda base (conda activate base)

ONLY FOR HAPLOIDS: USE THIS SCRIPT TO CHANGE CHROMOPAINTER INPUT TO HAPLOID 

for i in { A B C D E F G H };do (awk 'NR == 1  { print $1 /2 }' LACI0${i}.chromopainter; sed '2,3!d' LACI0${i}.chromopainter;sed '1,3d' LACI0${i}.chromopainter|sed  '0~2d')| cat > LACI0${i}.chromopainter.haploid;done


#Create recombination file for each chromosome (edit line 47 in the perl script makeuniformrecfile.pl to give a constant value of 4/1,000,000 (0.000004)) which means that in all the following calculations we will use a constant recombination rate of 0.4 cM/Kb (which is the average in S. cerevisiae (Cubillos et al, 2011)). 
Iterate this perl script over the files created in step 5 
Obs: makeuniformrecfile is inside finestructure folder

for i in A B C D E F G H;do perl ./fs_4.1.1/makeuniformrecfile.pl  LACI0${i}.chromopainter.haploid  LACI0${i}_rec.chromopainter.haploid;done

#Run chromopainter (v2) in mode ALLvsALL

for i in A B C D E F G H;do ./fs_4.1.1/fs_linux_glibc2.3 chromopainter -g LACI0${i}.chromopainter.haploid -r LACI0${i}_rec.chromopainter.haploid -t idfile.txt -o cp_LACI0${i} -a 0 0 -j 1;done 

#Run chromocombine to combine all chromopainter per chromosome outputs into a combined genome output (put all of them in a folder e.g. datos)

mkdir datos
mv cp_* datos

../fs_4.1.1/fs_linux_glibc2.3 chromocombine -d datos

#Run finestructure on the combined dataset

./fs_4.1.1/fs_linux_glibc2.3 finestructure  -x 100000 -y 100000 -z 1000 output.chunkcounts.out out.cidri.mcmc.xml


./fs_4.0.1/fs_linux_glibc2.3 finestructure  -x 100000 -k 2 -m T -t 1000000 output.chunkcounts.out structure_result.xml structure_tree.out
./fs_4.0.1/fs_linux_glibc2.3 finestructure  -e meancoincidence  output.chunkcounts.out structure_result.xml structure_meancoincidence.csv
./fs_4.0.1/fs_linux_glibc2.3 finestructure  -e X2  output.chunkcounts.out structure_result.xml structure_meanstate.csv



# --- fastStructure ---
# Variational Bayes approach to population structure (faster than STRUCTURE)
# Run via Docker to avoid dependency conflicts

sudo  docker run -it --rm -v $PWD:/data:rw  --platform=linux/amd64  fischuu/faststructure  bash 

for k in {2..10}; do python /fastStructure-1.0/structure.py -K $k --i=ldfilter --o=test_k$k.log --cv=1 ; done

# =============================================================================
# PART 4: GENE FLOW ANALYSIS — TreeMix
# =============================================================================
# TreeMix builds a population graph with migration edges from allele frequency data



vcftools --vcf filmax_all.recode.vcf --keep keep.txt --non-ref-ac-any 1 --recode-INFO-all --recode --out treemix_vcf_f1

plink2 --vcf treemix_vcf_f1.recode.vcf --double-id --allow-extra-chr --set-all-var-ids @:# --indep-pairwise 50 10 0.5 --out treemix_ldfilter
sed 's/polished:/polished\t/' treemix_ldfilter.prune.in > treemix_ldfilter.prune.in.vcftools
vcftools --vcf treemix_vcf_f1.recode.vcf --positions treemix_ldfilter.prune.in.vcftools --recode-INFO-all --recode --out treemix_vcf_f2

##Use STACKS to create a TREEMIX file from the filtered VCF. Here we use the popstacks.txt file we created before
populations -V treemix_vcf_f2.recode.vcf -O treemix -M popstacks.txt --treemix
##Delete the first line of the STACKS output, and gzip it
sed -i '1d' treemix_vcf_f2.recode.p.treemix
gzip treemix_vcf_f2.recode.p.treemix


# =============================================================================
# PART 5: GENE ORTHOLOGY — OrthoFinder
# =============================================================================
# Identify orthologs across Nanopore-assembled genomes using protein sequences
# Generated from GFF3 annotations with gffread

# Generate .faa protein FASTA files from assembly + GFF3 annotation


for SAMPLE in ${cat cepas.txt};do gffread ./ONT_seq/${SAMPLE}.final.gff3 -g ./ONT_seq/${SAMPLE}.assembly.final.fa -y ./protein2/${SAMPLE}.faa ; done


gffread ./ONT_seq/CDFM21L.1.final.gff3 -g ./ONT_seq/CDFM21L.1.assembly.final.fa -y ./protein/CDFM21L.faa

# Run OrthoFinder via Docker (all-vs-all DIAMOND search + orthogroup inference)
docker run -it --rm -v $PWD:/data:rw  --platform linux/amd64 staphb/orthofinder bash

# Inside Docker, run OrthoFinder on all .faa files (10 threads)
orthofinder -f . -t 10

# Remove mitochondrial chromosomes from GFF3 files before generating .faa
# (mitochondrial genes can cause issues with ortholog assignment)
grep -v "chrMT" ./ONT_seq/CDFM21L.1.final.gff3 > ./ONT_seq/CDFM21L.1.filtered.gff3

grep -v "chrMT" ./ONT_seq/aa4.final.gff3 > ./ONT_seq/aa4.filtered.gff3
gffread ./ONT_seq/aa4.filtered.gff3 -g ./ONT_seq/aa4.assembly.final.fa -y ./protein2/aa4.faa


grep -v "chrMT" ./ONT_seq/CL248.final.gff3 > ./ONT_seq/CL248.filtered.gff3
gffread ./ONT_seq/CL248.filtered.gff3 -g ./ONT_seq/CL248.assembly.final.fa -y ./protein2/CL248.faa






grep -v -E "chrMT|contig_26_1_0-42254_0_1_0-42302_0_1_0-42296_0|contig_12_1_0-74017_0_1_0-74100_0_1_0-74110_0|contig_18_1_0-14205_0_1_0-14506_0_1_0-14228_0|contig_23_1_0-14014_0_1_0-14024_0_1_0-14022_0" ./ONT_seq/QC18.final.gff3 > ./ONT_seq/QC18.filtered.gff3









